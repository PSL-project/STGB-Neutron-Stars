import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from stgb.figures import segment

def test_missing_stars_are_not_connected():
    fig,ax=plt.subplots()
    rows=[{'x':0,'y':1},{'x':1,'y':float('nan')},{'x':2,'y':3},{'x':3,'y':4}]
    assert segment(ax,rows,'x','y','tested',color='k')
    assert len(ax.lines)==2
    assert len(ax.lines[0].get_xdata())==1
    assert len(ax.lines[1].get_xdata())==2
    plt.close(fig)
