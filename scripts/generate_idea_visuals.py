"""Reproduce the explanatory figures and GIFs from computed road geometry.

Requires numpy, scipy, matplotlib, shapely, Pillow. Runs without the repository's
OCP dependencies. See README.md for provenance, scope and numerical checks.
"""
from pathlib import Path
import io
import os
import json
import argparse
os.environ.setdefault('MPLCONFIGDIR', '/tmp/ocp_fov_presentation_mpl')
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PatchPolygon, Patch
from shapely.geometry import Point, Polygon
from PIL import Image
from geometry_model import R, H, B, wd, cross, events

OUT=Path(__file__).resolve().parent
ANIM_OUT=OUT
BG='#F8FAFC'; INK='#17324D'; MUTED='#62768B'; GRID='#DFE7EF'
ROAD='#E5ECF2'; VISIBLE='#B7DDCE'; BLIND='#ED8C7F'; BLUE='#1768AB'
TEAL='#007E83'; GOLD='#CF7B10'; PURPLE='#8763AB'; BORDER='#3D5268'
COL=[GOLD,BLUE,PURPLE]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
    'text.color':INK,'axes.labelcolor':MUTED,'xtick.color':MUTED,
    'ytick.color':MUTED,'axes.edgecolor':GRID,'axes.spines.top':False,
    'axes.spines.right':False,'figure.facecolor':BG,'axes.facecolor':BG,
    'svg.fonttype':'none','savefig.facecolor':BG})
EV=events()
assert len(EV)==3
S=np.linspace(0,H,6001)
ROAD_POLY=Polygon(np.vstack([B(S,'L'),B(S[::-1],'R')]))
POCKET=Polygon(B(np.linspace(EV[0]['st'],EV[0]['first'][0],3001),'L'))
BEST=min((e for e in EV if e['side']!=e['first'][1]),key=lambda e:e['first'][0])
F=BEST['first'][0]
NEAR=Polygon(np.vstack([B(np.linspace(0,BEST['st'],3001),'L'),
                       B(np.linspace(F,0,4001),'R')]))
VISIBLE_POLY=NEAR.difference(POCKET)
FAR=ROAD_POLY.difference(NEAR)

def poly(ax,geom,fill,alpha=1,z=1,hatch=None):
    if geom.is_empty:return
    if geom.geom_type=='Polygon':
        ax.add_patch(PatchPolygon(np.array(geom.exterior.coords),facecolor=fill,
                                 edgecolor='#AF4C40' if hatch else 'none',alpha=alpha,
                                 hatch=hatch,linewidth=.3,zorder=z))
    elif hasattr(geom,'geoms'):
        for g in geom.geoms:poly(ax,g,fill,alpha,z,hatch)

def xylabel(ax,point,text,offset=(8,8),color=INK,fs=11):
    ax.annotate(text,point,xytext=offset,textcoords='offset points',fontsize=fs,
        fontweight='bold',color=color,zorder=20,
        bbox=dict(boxstyle='round,pad=.2',facecolor=BG,edgecolor='none',alpha=.86))

def base_road(ax,visible=False,pocket=False,stations=False):
    poly(ax,ROAD_POLY,ROAD)
    if visible:poly(ax,NEAR,VISIBLE)
    if pocket:poly(ax,POCKET,BLIND)
    left=B(S,'L');right=B(S,'R')
    ax.plot(left[:,0],left[:,1],color=BORDER,lw=1.6,zorder=4)
    ax.plot(right[:,0],right[:,1],color=BORDER,lw=1.6,zorder=4)
    for s in [0,H]:
        v=np.vstack([B(s,'L'),B(s,'R')]);ax.plot(v[:,0],v[:,1],color=BORDER,lw=1.1)
    center=np.column_stack([R*np.sin(S/R),R-R*np.cos(S/R)])
    ax.plot(center[:,0],center[:,1],color=MUTED,lw=.85,ls=(0,(4,5)),alpha=.7)
    if stations:
        for s in [2,6,10,14,18]:
            v=np.vstack([B(s,'L'),B(s,'R')]);ax.plot(v[:,0],v[:,1],color=MUTED,lw=.65,alpha=.4)
    ax.scatter([0],[0],s=140,marker='*',color=INK,zorder=20)
    xylabel(ax,np.zeros(2),'P',(-21,-17),fs=12)
    ax.set_xlim(-.8,12.8);ax.set_ylim(-1.4,14.)
    ax.set_aspect('equal');ax.set_xlabel('x [m]');ax.set_ylabel('y [m]')
    ax.set_xticks([0,4,8,12]);ax.set_yticks([0,4,8,12]);ax.grid(alpha=.35,lw=.55,zorder=0)

def event_cut(ax,j,selected=False,all_hit=False,alpha=1,labels=True):
    e=EV[j];T=e['T'];su,side=e['first'];U=B(su,side);color=COL[j]
    ax.plot([0,T[0]],[0,T[1]],color=color,ls=(0,(3,4)),lw=1.2,alpha=alpha,zorder=5)
    ax.plot([T[0],U[0]],[T[1],U[1]],color=color,lw=3.2 if selected else 2,
            ls='-' if j<2 else (0,(3,2)),alpha=alpha,zorder=6)
    if all_hit and j==0:
        Uo=B(e['hits']['R'][0],'R')
        ax.plot([U[0],Uo[0]],[U[1],Uo[1]],color=color,lw=1,ls=':',alpha=.65)
        ax.scatter(Uo[0],Uo[1],s=55,marker='o',facecolor=BG,edgecolor=color,zorder=10)
        if labels:xylabel(ax,Uo,r'$U_{1,\mathrm{opp}}$',(8,0),color,10)
    ax.scatter(T[0],T[1],s=65,color=color,edgecolor=BG,linewidth=.8,zorder=10)
    ax.scatter(U[0],U[1],s=60,marker='s',color=color,edgecolor=BG,linewidth=.8,zorder=10)
    if labels:
        offsets=[(-30,4),(-26,-24),(-34,9)]
        xylabel(ax,T,f'$T_{j+1}$',offsets[j],color)
        xylabel(ax,U,f'$U_{j+1}$',[(3,12),(8,-7),(9,3)][j],color)

def header(fig,kicker,title,subtitle):
    fig.text(.045,.969,kicker,fontsize=10,fontweight='bold',color=BLUE,va='top')
    fig.text(.045,.928,title,fontsize=21,fontweight='bold',va='top')
    fig.text(.045,.878,subtitle,fontsize=10.8,color=MUTED,va='top')
    fig.text(.045,.022,'ILLUSTRATIVE GEOMETRY  |  Smooth embedded road strip; computed events and ray intersections.',
             fontsize=8,color=MUTED)

def save_static(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=240)
    svg_path = OUT/(name+'.svg')
    fig.savefig(svg_path)
    # Keep generated SVG paths free of trailing whitespace in source control.
    svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines())+'\n')
    plt.close(fig)

def angle_figure(scan=None):
    fig=plt.figure(figsize=(14,8.4))
    header(fig,'01  /  EXTRACT THE EVENTS','Local angular extrema reveal the critical rays',
           'Both boundaries use the same road station s. Only forward, inward tangencies become cut candidates.')
    ax=fig.add_axes([.055,.125,.42,.69]);base_road(ax)
    ar=fig.add_axes([.56,.33,.385,.43])
    ss=np.linspace(.02,H,5000)
    for side,col in [('L',BLUE),('R',TEAL)]:
        q=B(ss,side);theta=np.degrees(np.arctan2(q[:,1],q[:,0]))
        ar.plot(ss,theta,color=col,lw=2,label=f'{side} boundary')
    ar.set_xlim(0,H);ar.set_ylim(-88,88)
    ar.set_xticks([0,3,6,9,12,15,18]);ar.set_yticks([-60,-30,0,30,60])
    ar.set_xlabel('Common road station s [m]');ar.set_ylabel('Bearing angle [deg]')
    ar.grid(color=GRID,lw=.65);ar.legend(loc='lower right',frameon=False,fontsize=10)
    for j,e in enumerate(EV):
        if scan is not None and e['st']>scan:continue
        T=e['T'];theta=np.degrees(np.arctan2(T[1],T[0]))
        ax.plot([0,T[0]],[0,T[1]],color=COL[j],lw=1.3,ls='--')
        ax.scatter(*T,s=68,color=COL[j],zorder=9)
        xylabel(ax,T,f'$T_{j+1}$',[(-28,8),(-27,-21),(-28,11)][j],COL[j])
        ar.scatter(e['st'],theta,s=70,color=COL[j],zorder=9)
        ar.annotate(f'$T_{j+1}$',(e['st'],theta),xytext=(0,16),textcoords='offset points',
                    color=COL[j],ha='center',fontweight='bold')
    if scan is not None:
        ar.axvline(scan,color=MUTED,lw=1,ls=':')
        for side,col in [('L',BLUE),('R',TEAL)]:
            q=B(scan,side);theta=np.degrees(np.arctan2(q[1],q[0]))
            ax.plot([0,q[0]],[0,q[1]],color=col,lw=1.5,alpha=.65)
            ax.scatter(*q,s=35,color=col,zorder=10)
            ar.scatter(scan,theta,s=45,color=col,zorder=10)
        stat=f'Scanning both boundaries at s = {scan:4.1f} m'
    else:stat='Three inward tangencies: T1, T2, T3'
    fig.text(.56,.247,stat,fontsize=12,fontweight='bold')
    fig.text(.56,.188,r'$\theta_b^{\prime}(s)=\det(B_b-P,\,B_b^{\prime})\,/\,\|B_b-P\|^2$',fontsize=16)
    fig.text(.56,.117,'At a nondegenerate extremum, the line of sight is tangent.\nThe other extrema do not pass the inward-event test.',
                fontsize=10.5,color=MUTED,linespacing=1.6)
    return fig

def comparison_figure(stage=5):
    fig=plt.figure(figsize=(14,8.4))
    titles=['Start from the tangent events','Compare both sides on the same ray',
            'A cross-road hit creates a front candidate','Retain later candidates for the second comparison',
            'Select the earliest global front','Classify locally. Select the front globally.']
    subtitles=['A shared station coordinate makes the boundary hits comparable.',
        'T1 returns to its own boundary first: a local blind pocket, not the global front.',
        'T2 reaches the opposite boundary first: it separates the road ahead.',
        'T3 is a later cross-road candidate. Hidden candidates can be redundant cuts.',
        'Take the minimum hit station among opposite-side candidates only.',
        'The scalar front location does not encode the local blind pocket.']
    header(fig,'02  /  TWO LEVELS OF COMPARISON',titles[stage],subtitles[stage])
    ax=fig.add_axes([.045,.115,.43,.7]);base_road(ax,visible=stage>=4,pocket=stage>=5)
    if stage>=4:poly(ax,FAR,'#CED7E1',alpha=.65,z=2)
    active=0 if stage==0 else min(stage,3)
    for j in range(active):event_cut(ax,j,selected=stage>=4 and j==1,all_hit=stage==1)
    if stage==0:
        for j,e in enumerate(EV):
            ax.scatter(*e['T'],s=65,color=COL[j],zorder=8)
            xylabel(ax,e['T'],f'$T_{j+1}$',(6,5),COL[j])
    if stage>=4:xylabel(ax,B(F,'R'),'FPV',(10,-23),BLUE,11)
    if stage>=5:
        ax.legend(handles=[Patch(facecolor=VISIBLE,label='Visible'),
                           Patch(facecolor=BLIND,label='Local blind pocket'),
                           Patch(facecolor='#D1DAE3',label='Beyond first front')],
                  loc='upper left',frameon=False,fontsize=9.5)
    fig.text(.53,.774,'1   CLASSIFY EACH EVENT',fontsize=12,fontweight='bold')
    tab=fig.add_axes([.52,.403,.45,.32]);tab.axis('off')
    tab.set_xlim(0,1);tab.set_ylim(0,1)
    cols=[.02,.18,.37,.57,.76]
    for x,head in zip(cols,['Event','s(T)','Same','Opposite','First hit']):
        tab.text(x,.93,head,fontsize=10,fontweight='bold')
    for j,e in enumerate(EV):
        yy=.7-j*.245
        tab.axhline(yy+.12,color=GRID,lw=.8)
        values=[f'T{j+1}',f"{e['st']:.3f}",
                f"{e['hits'][e['side']][0]:.3f}" if e['hits'][e['side']] else r'$\infty$',
                f"{e['hits']['R'][0]:.3f}",
                'Local pocket' if j==0 else 'Global front']
        for x,v in zip(cols,values):
            tab.text(x,yy,v,color=COL[j] if j<active else '#AFBBC7',
                     fontsize=10.8,fontweight='bold' if x==cols[0] else 'normal')
    fig.text(.53,.384,'All values are common road station s [m].\n'+r'$\infty$: no same-side return within the displayed window.',
             fontsize=9.2,color=MUTED,linespacing=1.5)
    fig.text(.53,.29,'2   SELECT THE GLOBAL FRONT',fontsize=12,fontweight='bold',
             color=INK if stage>=4 else '#AFBBC7')
    fig.text(.53,.227,r'$F=\min(10.306,\ 14.955)=10.306\ \mathrm{m}$',
             fontsize=16.5,color=BLUE if stage>=4 else '#AFBBC7')
    fig.text(.53,.16,'T1 contributes a pocket. T2 determines the global front.\nT3 is a later, redundant front candidate.',
             fontsize=11,color=INK if stage>=5 else MUTED,linespacing=1.6)
    fig.text(.53,.079,'Two comparison levels, not two scalar operations.',fontsize=9.5,color=MUTED)
    return fig

def partition_figure():
    fig=plt.figure(figsize=(13.8,8.5))
    header(fig,'03  /  COMPLETE VISIBILITY','The front is only one part of the answer',
           'A local blind pocket can lie before the maximum visible station and must remain excluded.')
    axs=[fig.add_axes([.055,.16,.37,.60]),fig.add_axes([.565,.16,.37,.60])]
    for k,ax in enumerate(axs):
        base_road(ax,visible=True,pocket=bool(k));poly(ax,FAR,'#CED7E1',.65,2)
        event_cut(ax,1,selected=True,labels=True)
        if k:
            event_cut(ax,0,labels=True)
            pc=np.array(POCKET.representative_point().coords)[0]
            ax.annotate('Local blind pocket',pc,xytext=(-.3,7.1),fontsize=10.5,
                color='#AF4C40',arrowprops=dict(arrowstyle='->',color='#AF4C40',lw=1.3),
                bbox=dict(boxstyle='round,pad=.3',facecolor=BG,edgecolor='none'),zorder=20)
            ax.legend(handles=[Patch(facecolor=VISIBLE,label='Visible'),
                               Patch(facecolor=BLIND,label='Local blind pocket'),
                               Patch(facecolor='#D1DAE3',label='Beyond first front')],
                      loc='upper left',frameon=False,fontsize=9)
        else:
            poly(ax,POCKET,'none',z=8,hatch='////')
            # Outline the pocket for the ablation without assigning it a visible label.
            a=np.array(POCKET.exterior.coords);ax.plot(a[:,0],a[:,1],color='#AF4C40',lw=1.4,zorder=9)
        ax.text(.02,1.065,['Front constraint alone','Full event partition'][k],transform=ax.transAxes,
                fontsize=14,fontweight='bold')
    fig.text(.064,.084,'The hatched pocket is not excluded by the front alone.',fontsize=10.5,color=MUTED)
    fig.text(.575,.084,'Visible = near-side region minus the local blind pocket.',fontsize=10.5,color=MUTED)
    fig.text(.49,.43,'→',fontsize=27,color=MUTED,ha='center')
    return fig

def image_from_fig(fig):
    buf=io.BytesIO();fig.savefig(buf,format='png',dpi=105)
    buf.seek(0);im=Image.open(buf).convert('RGB').copy();plt.close(fig);return im

def write_gif(frames,name,durations):
    frames[0].save(ANIM_OUT/name,save_all=True,append_images=frames[1:],duration=durations,
                   loop=0,optimize=False,disposal=2)

def exact_margin(q):
    # Continuous ray oracle, independent of the event/crosscut construction.
    def clearance(lam):
        lam=np.asarray(lam);x=lam*q[0];y=lam*q[1]
        s=R*np.arctan2(x,R-y);r=np.hypot(x,R-y)
        return np.minimum(r-(R-wd(s,'L')[0]),(R+wd(s,'R')[0])-r)
    ls=np.linspace(0,1,301);ms=clearance(ls)
    candidate=[float(ms[0]),float(ms[-1])]
    for i in np.flatnonzero((ms[1:-1]<=ms[:-2])&(ms[1:-1]<=ms[2:]))+1:
        opt=minimize_scalar(clearance,bounds=(ls[i-1],ls[i+1]),method='bounded',
                            options={'xatol':1e-13})
        candidate.append(float(opt.fun))
    return min(candidate)

def verify():
    rng=np.random.default_rng(20260924);failures=[];skipped=0;count=4000
    for _ in range(count):
        s=rng.uniform(.0001,H);wl=wd(s,'L')[0];wr=wd(s,'R')[0]
        n=rng.uniform(-wr,wl);r=R-n
        q=np.array([r*np.sin(s/R),R-r*np.cos(s/R)])
        margin=exact_margin(q)
        if abs(margin)<3e-5:skipped+=1;continue
        pred=VISIBLE_POLY.covers(Point(q))
        if pred!=(margin>=0):failures.append({'s':s,'n':n,'margin':margin,'predicted':pred})
    details=[]
    for j,e in enumerate(EV):
        su,side=e['first'];U=B(su,side);T=e['T']
        ls=np.linspace(.002,.999,1001);q=T[None,:]+ls[:,None]*(U-T)
        station=R*np.arctan2(q[:,0],R-q[:,1]);r=np.hypot(q[:,0],R-q[:,1])
        margin=np.minimum(r-(R-wd(station,'L')[0]),R+wd(station,'R')[0]-r)
        details.append({'event':f'T{j+1}','tangent_station':e['st'],'same_hit':e['hits']['L'][0] if e['hits']['L'] else None,
            'opposite_hit':e['hits']['R'][0] if e['hits']['R'] else None,
            'classification':'same-side pocket' if side=='L' else 'opposite-side front',
            'tangency_residual':float(abs(cross(T,B(e['st'],'L',True)))),
            'hit_alignment_residual':float(abs(cross(U,T))),
            'cut_minimum_sampled_interior_margin':float(margin.min()),
            'cut_station_strictly_increasing':bool(np.all(np.diff(station)>0))})
        assert margin.min()>0 and np.all(np.diff(station)>0)
    result={'scope':'Numerical checks of an illustrative smooth road, not verification of the original source pipeline.',
        'road':{'R':R,'H':H,'left_width_range':[float(wd(S,'L')[0].min()),float(wd(S,'L')[0].max())],
                'right_width_range':[float(wd(S,'R')[0].min()),float(wd(S,'R')[0].max())],
                'injectivity_reason':'polar angle s/R in [0,1.8] and all radii positive; each station is a distinct radial fibre',
                'polygon_is_valid':bool(ROAD_POLY.is_valid)},
        'events':details,'selected_FPV_station':F,'visible_area_m2':float(VISIBLE_POLY.area),
        'local_pocket_area_m2':float(POCKET.area),'independent_oracle':{
            'method':'Continuous line-of-sight radial clearances; dense bracketing plus local bounded minimization',
            'seed':20260924,'samples':count,'near_boundary_skipped':skipped,
            'mismatches':len(failures),'examples':failures[:5]},
        'limitations':['Numerical event enumeration and oracle checks are not a formal certificate.',
                       'No actual-track performance, OCP trajectory, sensor noise or speed improvement is claimed.',
                       'Figures assume opaque road exterior and no internal occluders.',
                       'T3 is intentionally retained as a hidden, redundant global-front candidate.']}
    (OUT/'geometry_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    assert ROAD_POLY.is_valid and not failures,result
    return result

def main():
    global OUT,ANIM_OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=OUT,help='Directory for PNG, SVG and geometry_checks.json')
    parser.add_argument('--animations',type=Path,help='Directory for GIFs; defaults to --out')
    args=parser.parse_args()
    OUT=args.out.resolve();ANIM_OUT=(args.animations or OUT).resolve()
    OUT.mkdir(parents=True,exist_ok=True);ANIM_OUT.mkdir(parents=True,exist_ok=True)
    result=verify()
    save_static(angle_figure(),'01_tangent_events')
    save_static(comparison_figure(),'02_two_stage_selection')
    save_static(partition_figure(),'03_complete_visibility')
    scanvalues=np.r_[np.linspace(.07,H,64),np.repeat(H,12)]
    scanframes=[image_from_fig(angle_figure(float(s))) for s in scanvalues]
    write_gif(scanframes,'01_boundary_angle_scan.gif',[140]*64+[300]*12)
    pipeline=[image_from_fig(comparison_figure(k)) for k in range(6)]
    write_gif(pipeline,'02_two_stage_story.gif',[1500,2300,1900,1900,2600,3400])
    print(json.dumps({'events':len(EV),'FPV':F,'oracle':result['independent_oracle'],
                      'files':[p.name for p in OUT.glob('*') if p.suffix in ['.png','.svg','.gif']]},indent=2))

if __name__=='__main__':main()
