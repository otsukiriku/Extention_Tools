import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

def plot_scatter(xdata:list, ydata:list, _figname = "test", xlabel='x', ylabel='y'):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    #ax.set_ylim(0, 10)
    ax.set_xlabel("Prove Radius(Aungstrom)", size = 16,)
    ax.set_ylabel("Coverage(%)", size = 16,)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    #plt.figure(figsize=(4,3))
    ax.scatter(xdata,ydata)
    test = [[i] for i in xdata]
    X = np.array(test)
    Y = np.array(ydata)
    Ir = LinearRegression()
    Ir.fit(X,Y)
    plt.plot(X,Ir.predict(X), color="red")
    ax.text(0.95,0.05,"coefficient="+str(Ir.coef_), transform=ax.transAxes, horizontalalignment="right")
    fig.savefig(_figname+".png",dpi=250)

def simpleLR_withscatter(xdata:list, ydata:list, _figname = "test", xlabel='x', ylabel='y'):
    """datax,y should be 2dimension list"""
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    #ax.set_ylim(0, 10)
    ax.set_xlabel("Prove Radius(Aungstrom)", size = 16,)
    ax.set_ylabel("Coverage(%)", size = 16,)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_xlim(5,55)
    ax.set_ylim(0,80)
    #plt.figure(figsize=(4,3))
    print("xdata dimension")
    print(len(xdata))
    if len(xdata) == 1:
        ax.scatter(xdata,ydata)
        test = [[i] for i in xdata]
        X = np.array(test)
        Y = np.array(ydata)
        Ir = LinearRegression()
        Ir.fit(X,Y)
        plt.plot(X,Ir.predict(X), color="red")
        ax.text(0.95,0.05,"coefficient="+str(Ir.coef_), transform=ax.transAxes, horizontalalignment="right")
        fig.savefig(_figname+".png",dpi=250)
        plt.show()
    else:
        for idx, (x,y) in enumerate(zip(xdata,ydata)):
            ax.scatter(x,y,c=dot_color[idx])
            test = [[i] for i in x]
            X = np.array(test)
            Y = np.array(y)
            Ir = LinearRegression()
            Ir.fit(X,Y)
            plt.plot(X,Ir.predict(X), color=line_color[idx])
            ax.text(0.95,0.05+0.05*idx, dot_color[idx]+":coefficient="+str(Ir.coef_), transform=ax.transAxes, horizontalalignment="right")
        fig.savefig(_figname+".png",dpi=250)


