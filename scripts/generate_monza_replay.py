#!/usr/bin/env python3
"""Render a full-lap Monza geometry replay, preserving the original GIFs.

The existing sampled event algorithm supplies the candidate and FPV markers.
An independent polygonal ray cast supplies the reference visibility shading.
Animation time is presentation time, not an optimized vehicle trajectory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import FancyBboxPatch, Polygon
import numpy as np
from PIL import Image

from monza_replay_geometry import prepare_replay

BG = '#080f1a'
PANEL = '#101c2b'
GRID = '#1c2b3d'
TEXT = '#edf5fd'
MUTED = '#8ca1b8'
CYAN = '#5ee5ef'
BLUE = '#338dca'
GOLD = '#ffca6b'
RED = '#ed7d96'
ROAD = '#2c3c50'


def local(points, p, heading):
    a = np.asarray(points, dtype=float).reshape(-1, 2) - p
    h = np.array([np.cos(heading), np.sin(heading)])
    right = np.array([h[1], -h[0]])
    return np.column_stack((a @ right, a @ h))


def set_line(artist, xy):
    xy = np.asarray(xy).reshape(-1, 2)
    artist.set_data(xy[:, 0], xy[:, 1])


class ReplayDashboard:
    def __init__(self, data, width=1600, height=1000):
        self.data = data
        self.frames = data['frames']
        self.length = float(data['length'])
        plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                             'text.color': TEXT, 'axes.labelcolor': MUTED,
                             'xtick.color': MUTED, 'ytick.color': MUTED,
                             'savefig.facecolor': BG, 'axes.unicode_minus': False})
        self.fig = plt.figure(figsize=(width / 100, height / 100), dpi=100, facecolor=BG)
        self.fig.text(.043, .95, 'MONZA', fontsize=31, weight='bold', color=TEXT)
        self.fig.text(.209, .95, '/  A LAP THROUGH VISIBILITY', fontsize=18, color=CYAN)
        self.fig.text(.044, .917, 'COMMON ROAD STATIONS  /  TANGENTIAL EVENTS  /  TWO-STAGE FRONT SELECTION',
                      fontsize=9.2, color=MUTED)
        self.fig.text(.956, .955, 'OCP FOV', fontsize=12, color=CYAN, ha='right', weight='bold')
        self.fig.text(.956, .923, 'FULL CIRCUIT  /  GEOMETRY REPLAY', fontsize=8.5, color=MUTED, ha='right')
        self.fig.add_artist(plt.Line2D([.04, .96], [.893, .893], transform=self.fig.transFigure,
                                      color=GRID, lw=1))
        self.fig.text(.05, .86, '01  /  CIRCUIT OVERVIEW', fontsize=10.2, color=MUTED, weight='bold')
        self.fig.text(.36, .86, '02  /  FOLLOW CAMERA', fontsize=10.2, color=MUTED, weight='bold')
        self.mode = self.fig.text(.953, .86, '', fontsize=10, color=GOLD, ha='right')
        for x, w in [(.035, .285), (.343, .622)]:
            self.fig.add_artist(FancyBboxPatch((x, .278), w, .555,
                boxstyle='round,pad=0.006,rounding_size=0.011',
                transform=self.fig.transFigure, facecolor=PANEL, edgecolor=GRID,
                linewidth=1, zorder=-1))

        self.map = self.fig.add_axes([.048, .30, .26, .51], facecolor=PANEL)
        center = data['whole_center']
        self.map.plot(center[:, 0], center[:, 1], color='#20364a', lw=9, solid_capstyle='round')
        self.map.plot(center[:, 0], center[:, 1], color='#6b8197', lw=2.2, solid_capstyle='round')
        self.map.plot(center[0, 0], center[0, 1], marker='s', color=TEXT, ms=4)
        self.map.annotate('START / FINISH', center[0], xytext=(12, -5),
                          textcoords='offset points', fontsize=7.5, color=MUTED)
        self.progress_glow, = self.map.plot([], [], color=CYAN, lw=8, alpha=.1, solid_capstyle='round')
        self.progress_line, = self.map.plot([], [], color=CYAN, lw=2.8, solid_capstyle='round')
        self.map_window, = self.map.plot([], [], color=GOLD, lw=4, alpha=.95, solid_capstyle='round')
        self.map_halo, = self.map.plot([], [], 'o', color=CYAN, ms=19, alpha=.13)
        self.map_car, = self.map.plot([], [], 'o', color=TEXT, ms=6, mec=CYAN, mew=1.5)
        self.map.set_aspect('equal')
        self.map.set_xlim(center[:, 0].min()-170, center[:, 0].max()+210)
        self.map.set_ylim(center[:, 1].min()-70, center[:, 1].max()+75)
        self.map.axis('off')
        self.fig.text(.053, .289, f'{self.length / 1000:.3f} km  /  PERIODIC ROAD', fontsize=9.3, color=MUTED)

        self.follow = self.fig.add_axes([.354, .30, .600, .51], facecolor=PANEL)
        ax = self.follow
        ax.set_xlim(-241, 241)
        ax.set_ylim(-49, 207)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xticks(np.arange(-200, 201, 50))
        ax.set_yticks(np.arange(-25, 201, 25))
        ax.grid(color=GRID, lw=.55, alpha=.8)
        ax.tick_params(length=0, labelbottom=False, labelleft=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        # A faint heading guide echoes the original animation; it is not a sensor model.
        ang = np.linspace(-.58, .58, 70)
        fan = np.vstack(([0, 0], np.column_stack((120*np.sin(ang), 120*np.cos(ang))), [0, 0]))
        ax.add_patch(Polygon(fan, facecolor=BLUE, alpha=.045, edgecolor='none', zorder=0))
        ax.plot(fan[1:-1, 0], fan[1:-1, 1], color=BLUE, alpha=.20, lw=.7)
        self.road = Polygon(np.zeros((3, 2)), facecolor=ROAD, edgecolor='none', zorder=1)
        self.window = Polygon(np.zeros((3, 2)), facecolor=RED, alpha=.36, edgecolor='none', zorder=2)
        self.visible = Polygon(np.zeros((3, 2)), facecolor=CYAN, alpha=.68, edgecolor='none', zorder=3)
        for patch in [self.road, self.window, self.visible]:
            ax.add_patch(patch)
        self.bound_l, = ax.plot([], [], color='#93a6bb', lw=1.25, zorder=4)
        self.bound_r, = ax.plot([], [], color='#93a6bb', lw=1.25, zorder=4)
        self.center, = ax.plot([], [], color='#bccddc', lw=.6, ls=(0, (6, 7)), alpha=.45, zorder=4)
        self.trail_glow, = ax.plot([], [], color=CYAN, lw=7, alpha=.11, zorder=4)
        self.trail, = ax.plot([], [], color=CYAN, lw=1.5, alpha=.75, zorder=5)
        self.rays = LineCollection([], colors=GOLD, linewidths=.8, alpha=.44, zorder=5)
        ax.add_collection(self.rays)
        self.candidates, = ax.plot([], [], 'o', ms=4.5, mfc=PANEL, mec=GOLD, mew=1.25, ls='none', zorder=6)
        self.front_glow, = ax.plot([], [], color=GOLD, lw=10, alpha=.16, zorder=6)
        self.front, = ax.plot([], [], color=GOLD, lw=2.9, zorder=7, solid_capstyle='round')
        self.front_points, = ax.plot([], [], 'o', mfc=GOLD, mec=PANEL, mew=1, ms=7, ls='none', zorder=8)
        self.itp_text = ax.text(0, 0, 'ITP', color=GOLD, fontsize=8.8, weight='bold', zorder=9,
                                bbox=dict(boxstyle='round,pad=.3', facecolor=PANEL, edgecolor='none', alpha=.94))
        self.fpv_text = ax.text(0, 0, 'FPV', color=GOLD, fontsize=8.8, weight='bold', zorder=9,
                                bbox=dict(boxstyle='round,pad=.3', facecolor=PANEL, edgecolor='none', alpha=.94))
        # Enlarged car symbol for readability; no vehicle dimensions are inferred.
        ax.scatter([0], [0], s=640, color=CYAN, alpha=.10, edgecolors='none', zorder=9)
        car = np.array([[-2.2, -4.7], [-2.2, 3.4], [-1.6, 5.6], [1.6, 5.6], [2.2, 3.4], [2.2, -4.7]])
        ax.add_patch(Polygon(car, facecolor=TEXT, edgecolor=CYAN, lw=1.1, zorder=10))
        ax.add_patch(Polygon([[-1.6,.8],[-1.3,3],[1.3,3],[1.6,.8]], facecolor=BLUE, edgecolor='none', zorder=11))
        ax.plot([-220, -170], [-34, -34], color=MUTED, lw=1.2)
        ax.plot([-220,-220],[-36,-32],color=MUTED,lw=1)
        ax.plot([-170,-170],[-36,-32],color=MUTED,lw=1)
        ax.text(-195,-44,'50 m',ha='center',fontsize=9.1,color=MUTED)
        ax.text(220,-37,'HEADING UP',ha='right',fontsize=9.1,color=MUTED)

        for x, color, label in [(.361, CYAN, 'VISIBLE REFERENCE'), (.526, RED, 'OCCLUDED REFERENCE'),
                                (.718, GOLD, 'SAMPLED EVENTS / FRONT')]:
            self.fig.add_artist(plt.Line2D([x,x+.017],[.818,.818],transform=self.fig.transFigure,color=color,lw=3))
            self.fig.text(x+.023,.813,label,fontsize=9.2,color=MUTED)

        self.stats = []
        for x, title in [(.048,'ROAD STATION'),(.298,'SAMPLED FPV AHEAD'),(.565,'DETECTED CANDIDATES'),(.808,'LAP PROGRESS')]:
            self.fig.text(x,.245,title,fontsize=8.3,color=MUTED,weight='bold')
            self.stats.append(self.fig.text(x,.212,'',fontsize=21,color=TEXT,weight='bold'))
        self.fig.text(.049,.168,'03  /  FRONT DISTANCE AROUND THE LAP',fontsize=8.5,color=MUTED,weight='bold')
        self.fig.text(.952,.168,'GAPS = NO FRONT DETECTED IN THE 200 m WINDOW',fontsize=9,color=MUTED,ha='right')
        self.chart = self.fig.add_axes([.06,.073,.885,.079],facecolor=BG)
        station = np.array([f['s0'] for f in self.frames])
        fpv = np.array([f['fpv_delta'] if f['fpv_delta'] is not None else np.nan for f in self.frames])
        self.chart.plot(station/1000,fpv,color=GOLD,lw=1.15,marker='.',ms=2.2,alpha=.85)
        self.chart.axhline(200,color=GRID,lw=.7,ls='--')
        self.chart.set_xlim(0,self.length/1000)
        self.chart.set_ylim(0,218)
        self.chart.set_yticks([0,100,200])
        self.chart.set_yticklabels(['0','100','200 m'],fontsize=9.2)
        self.chart.set_xticks(np.arange(0,6))
        self.chart.set_xticklabels([f'{v} km' for v in range(6)],fontsize=9.2)
        self.chart.tick_params(length=0,pad=6)
        for spine in self.chart.spines.values():
            spine.set_visible(False)
        self.cursor = self.chart.axvline(0,color=CYAN,lw=1.25,alpha=.8)
        self.chart_dot, = self.chart.plot([],[],'o',color=GOLD,ms=5)
        self.fig.text(.045,.019,'MONZA.CSV  /  200 m WINDOW  /  SAMPLED GEOMETRY',
                      fontsize=9.2,color=MUTED)
        self.fig.text(.956,.019,'GEOMETRY REPLAY  /  NOT A SPEED TRACE',fontsize=9.2,color=MUTED,ha='right')

    def draw(self, index):
        f = self.frames[index]
        p, heading = np.asarray(f['p']), f['heading']
        world_history = np.asarray([e['p'] for e in self.frames[:index+1]])
        for artist in [self.progress_line,self.progress_glow]:
            set_line(artist,world_history)
        set_line(self.map_window,f['forward_center'])
        for artist in [self.map_car,self.map_halo]:
            artist.set_data([p[0]],[p[1]])
        cl, cr = local(f['context_left'],p,heading), local(f['context_right'],p,heading)
        self.road.set_xy(np.vstack((cl,cr[::-1])))
        self.window.set_xy(local(f['window'],p,heading))
        self.visible.set_xy(local(f['visible'],p,heading))
        set_line(self.bound_l,cl)
        set_line(self.bound_r,cr)
        set_line(self.center,local(f['context_center'],p,heading))
        trail = local(world_history[max(0,index-14):],p,heading)
        for artist in [self.trail,self.trail_glow]:
            set_line(artist,trail)
        tps = local(f['tangent_points'],p,heading)
        set_line(self.candidates,tps)
        self.rays.set_segments([np.array([[0.,0.],pt]) for pt in tps])
        front = local(f['border'],p,heading)
        for artist in [self.front,self.front_glow,self.front_points]:
            set_line(artist,front)
        has_front = len(front)>=2
        for label in [self.itp_text,self.fpv_text]:
            label.set_visible(has_front)
        if has_front:
            self.itp_text.set_position(front[0]+[-16,9])
            self.fpv_text.set_position(front[-1]+[8,-9])
        mismatch = f.get('front_reference_agrees') is False
        self.front.set_linestyle((0, (3, 2)) if mismatch else '-')
        self.mode.set_text(('SAMPLED FRONT / REFERENCE DIFFERS' if mismatch else 'SAMPLED FRONT DETECTED') if has_front else 'NO FRONT DETECTED WITHIN WINDOW')
        self.mode.set_color(GOLD if has_front else MUTED)
        self.stats[0].set_text(f"{f['s0']/1000:05.2f} km")
        self.stats[1].set_text(f"{f['fpv_delta']:.1f} m" if has_front else 'NOT DETECTED')
        self.stats[1].set_fontsize(21 if has_front else 15)
        self.stats[1].set_color(GOLD if has_front else MUTED)
        self.stats[2].set_text(f'{len(tps):02d}')
        self.stats[3].set_text(f"{100*f['s0']/self.length:05.1f}%")
        self.cursor.set_xdata([f['s0']/1000]*2)
        if has_front:
            self.chart_dot.set_data([f['s0']/1000],[f['fpv_delta']])
        else:
            self.chart_dot.set_data([],[])
        self.fig.canvas.draw()
        return np.asarray(self.fig.canvas.buffer_rgba())[:,:,:3].copy()


def encoder_executable():
    import imageio_ffmpeg
    return imageio_ffmpeg.get_ffmpeg_exe()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--track',type=Path,default=Path('Monza.csv'))
    ap.add_argument('--out',type=Path,default=Path('docs/animations'))
    ap.add_argument('--figures',type=Path,default=Path('docs/figures'))
    ap.add_argument('--frames',type=int,default=720)
    ap.add_argument('--fps',type=int,default=24)
    ap.add_argument('--width',type=int,default=1600)
    ap.add_argument('--height',type=int,default=1000)
    ap.add_argument('--preview-only',action='store_true')
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    args.figures.mkdir(parents=True,exist_ok=True)
    data=prepare_replay(args.track,frames=args.frames,look=200.,step=.5)
    dashboard=ReplayDashboard(data,args.width,args.height)
    frames=data['frames']
    interesting=[i for i,f in enumerate(frames) if f['fpv_delta'] is not None and f.get('front_reference_agrees') is not False]
    poster=min(interesting,key=lambda i:abs(frames[i]['s0']-1045.45)) if interesting else len(frames)//5
    preview_indices=sorted(set([0,len(frames)//6,len(frames)//3,len(frames)//2,2*len(frames)//3,5*len(frames)//6,len(frames)-1,poster]))
    if args.preview_only:
        for i in preview_indices:
            frame=dashboard.draw(i)
            Image.fromarray(frame).save(args.figures/f'monza_preview_{i:04d}.png')
    Image.fromarray(dashboard.draw(poster)).save(args.figures/'monza_full_lap_poster.png')
    report=dict(data.get('audit',{}))
    report.update({'frames':len(frames),'fps':args.fps,'duration_seconds':len(frames)/args.fps,
                   'video_pixels':[args.width,args.height], 'poster_frame':poster,
                   'source_csv_sha256':hashlib.sha256(args.track.read_bytes()).hexdigest()})
    if args.preview_only:
        (args.figures/'monza_replay_checks.json').write_text(json.dumps(report,indent=2))
        print('Preview frames ready.',flush=True)
        return
    ffmpeg=encoder_executable()
    video=args.out/'monza_full_lap.mp4'
    proc=subprocess.Popen([ffmpeg,'-y','-loglevel','error','-f','rawvideo','-vcodec','rawvideo',
        '-s',f'{args.width}x{args.height}','-pix_fmt','rgb24','-r',str(args.fps),'-i','-',
        '-an','-c:v','libx264','-preset','medium','-crf','20','-pix_fmt','yuv420p',
        '-movflags','+faststart',str(video)],stdin=subprocess.PIPE)
    try:
        for i in range(len(frames)):
            proc.stdin.write(dashboard.draw(i).tobytes())
            if i%60==0:
                print(f'Rendered {i+1}/{len(frames)} frames',flush=True)
        proc.stdin.close()
        if proc.wait()!=0:
            raise RuntimeError('Video encoding failed')
    finally:
        if proc.poll() is None:
            proc.terminate()
    gif=args.out/'monza_full_lap.gif'
    # A single shared palette prevents color flicker in the looping README GIF.
    gif_fps=args.fps/2
    vf=f'fps={gif_fps},scale=960:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128:stats_mode=full[p];[b][p]paletteuse=dither=bayer:bayer_scale=4'
    subprocess.run([ffmpeg,'-y','-loglevel','error','-i',str(video),'-filter_complex',vf,'-loop','0',str(gif)],check=True)
    with Image.open(gif) as im:
        gif_frames=im.n_frames
        duration=0
        for i in range(gif_frames):
            im.seek(i);duration+=im.info.get('duration',0)
        report.update({'gif_frames':gif_frames,'gif_pixels':list(im.size),'gif_duration_seconds':duration/1000,
                       'gif_loop':im.info.get('loop',0)})
    report['outputs_bytes']={p.name:p.stat().st_size for p in [video,gif]}
    (args.figures/'monza_replay_checks.json').write_text(json.dumps(report,indent=2))
    plt.close(dashboard.fig)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
