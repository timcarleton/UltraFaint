from numpy import log10,shape,where,random,array,cumsum,repeat
def getvmaxfinal(mloss,gamma,minrmax=False,withscatter=False):
    if minrmax:
        if gamma>0.75:
            mu=0.12
            nu=0.29
        else:
            mu=0.53
            nu=0.53
    else:
        if gamma>1.25 and gamma<1.75:
            mu=0.4
            nu=0.24
        elif gamma>0.75:
            #using new fits from Raphael
            #mu=0.67
            #nu=0.37
            mu=0.4
            nu=0.3
        elif gamma>0.25:
            mu=0.4
            nu=0.35
        else:
            mu=0.4
            nu=0.37
    if withscatter and gamma<.25:
        mlbins=array([10,.5,.2,.1,0])
        histx=[array([1,0]),array([1,.85]),array([1,.8]),array([1,.8])]
        histy=[array([1,0]),array([.9,.1]),array([.625,.375]),array([.5,.5])]
        if len(shape(mloss))==0:
            w=where(mloss>mlbins)[0][0]-1
            r=random.rand()
            wx=where(r<cumsum(histy[w]))[0][0]
            return 2**mu*mloss**nu/(1+mloss)**mu*histx[w][wx]
        else:
            ret=[]
            for i in mloss:
                w=where(i>mlbins)[0][0]-1
                r=random.rand()
                wx=where(r<cumsum(histy[w]))[0][0]
                ret.append(2**mu*i**nu/(1+i)**mu*histx[w][wx])
            return array(ret)


    return 2**mu*mloss**nu/(1+mloss)**mu


def getrmaxfinal(mloss,gamma,minrmax=False,withscatter=False):
    if minrmax:
        if gamma>0.75:
            mu=-0.17
            nu=0.43
        else:
            mu=-0.85
            nu=0.023
    else:
        if gamma>1.25 and gamma<1.75:
            mu=0
            nu=0.48
        elif gamma>0.75:
            #using new fits from Raphael
            #mu=0.95
            #nu=0.58
            mu=-0.3
            nu=0.4
        elif gamma>0.25:
            mu=-0.4
            nu=0.27
        else:
            mu=-1.3
            nu=0.05

    if withscatter and gamma<.25:
        mlbins=array([10,.5,.2,.1,0])
        histx=[array([1,0]),array([1,1.1]),array([1,1.2]),array([1,1.2])]
        histy=[array([1,0]),array([.9,.1]),array([.625,.375]),array([.5,.5])]
        if len(shape(mloss))==0:
            w=where(mloss>mlbins)[0][0]-1
            r=random.rand()
            wx=where(r<cumsum(histy[w]))[0][0]
            return 2**mu*mloss**nu/(1+mloss)**mu*histx[w][wx]
        else:
            ret=[]
            for i in mloss:
                w=where(i>mlbins)[0][0]-1
                r=random.rand()
                wx=where(r<cumsum(histy[w]))[0][0]
                ret.append(2**mu*i**nu/(1+i)**mu*histx[w][wx])
            return array(ret)

    return 2**mu*mloss**nu/(1+mloss)**mu

def getrhfinal(mloss,gamma,rovera=0.1,minrmax=False,withscatter=False):


    if minrmax:
        
        if gamma>0.5:
            
            alpha2=0.075
            beta2=-0.037

            alpha1=-0.10
            beta1=-0.067

        else:

            alpha2=0.019
            beta2=-0.18

            alpha1=0.0036
            beta1=-0.20

    else:

        if gamma>0.5:
            
            alpha2=1.22
            beta2=0.33

            alpha1=1.49
            beta1=0.35

        else:

            alpha2=1.63
            beta2=0.03
            alpha1=2.91
            beta1=0.15

    r2=2**alpha2*mloss**beta2/(1+mloss)**alpha2
    r1=2**alpha1*mloss**beta1/(1+mloss)**alpha1

    interpvalue=10**(log10(r2)+(log10(r2)-log10(r1))/(log10(.2)-log10(.1))*(log10(rovera)-log10(.2)))


    if len(shape(rovera))==0 and len(shape(mloss))!=0:
        rovera=repeat(rovera,len(mloss))

    if withscatter and gamma<.25:
        mlbins=array([10,.5,.2,.1,0])
        histx1=[array([1,0]),array([1,1.3]),array([1,1.4]),array([1,1.4])]
        histx2=[array([1,0]),array([1,1.25]),array([1,1.5]),array([1,1.5])]
        histy=[array([1,0]),array([.9,.1]),array([.625,.375]),array([.5,.5])]
        if len(shape(mloss))==0:
            w=where(mloss>mlbins)[0][0]-1
            r=random.rand()
            wx=where(r<cumsum(histy[w]))[0][0]
            return interpvalue*10**(log10(histx2[w][wx])+(log10(histx2[w][wx])-log10(histx1[w][wx]))/(log10(.2)-log10(.1))*(log10(rovera)-log10(.2)))
        else:
            ret=[]
            for i in range(len(mloss)):
                w=where(mloss[i]>mlbins)[0][0]-1
                r=random.rand()
                wx=where(r<cumsum(histy[w]))[0][0]
                ret.append(interpvalue[i]*10**(log10(histx2[w][wx])+(log10(histx2[w][wx])-log10(histx1[w][wx]))/(log10(.2)-log10(.1))*(log10(rovera[i])-log10(.2))))
            return array(ret)

    return interpvalue

def getmstarfinal(mloss,gamma,rovera=0.1,minrmax=False,withscatter=False):

    if minrmax:

        if gamma>0.5:

            alpha2=0.53
            beta2=0.20

            alpha1=0.15
            beta1=0.051
        else:

            alpha2=0.33
            beta2=0.14
            alpha1=0.14
            beta1=0.053
    else:
    
        if gamma>0.5:

            alpha2=3.57
            beta2=2.06
            alpha1=3.43
            beta1=1.86

        else:

            alpha2=0.82
            beta2=0.82
            alpha1=1.43
            beta1=0.69

    m1=2**alpha1*mloss**beta1/(1+mloss)**alpha1
    m2=2**alpha2*mloss**beta2/(1+mloss)**alpha2

    interpvalue=10**(log10(m2)+(log10(m2)-log10(m1))/(log10(.2)-log10(.1))*(log10(rovera)-log10(.2)))

    if len(shape(rovera))==0 and len(shape(mloss))!=0:
        rovera=repeat(rovera,len(mloss))
    if withscatter and gamma<.25:
        mlbins=array([10,.5,.2,.1,0])
        histx1=[array([1,0]),array([1,.9]),array([1,.95]),array([1,.95])]
        histx2=[array([1,0]),array([1,.7]),array([1,.95]),array([1,.95])]
        histy=[array([1,0]),array([.9,.1]),array([.625,.375]),array([.5,.5])]
        if len(shape(mloss))==0:
            w=where(mloss>mlbins)[0][0]-1
            r=random.rand()
            wx=where(r<cumsum(histy[w]))[0][0]
            return interpvalue*10**(log10(histx2[w][wx])+(log10(histx2[w][wx])-log10(histx1[w][wx]))/(log10(.2)-log10(.1))*(log10(rovera)-log10(.2)))
        else:
            ret=[]
            for i in range(len(mloss)):
                w=where(mloss[i]>mlbins)[0][0]-1
                r=random.rand()
                wx=where(r<cumsum(histy[w]))[0][0]
                ret.append(interpvalue[i]*10**(log10(histx2[w][wx])+(log10(histx2[w][wx])-log10(histx1[w][wx]))/(log10(.2)-log10(.1))*(log10(rovera[i])-log10(.2))))
            return array(ret)

    
    return interpvalue

def getmstarfinalnew(mlossrmax,gamma,rovera=0.1):



    return 2**alpha*mlossrmax**beta/(1+mlossrmax)**alpha
    
def getmstarfinaltim(mlossrmax,gamma,rovera=0.1):

    rmaxratio=getrmaxfinal(mlossrmax,gamma,minrmax=True)
    return getmstarfinal((rovera+1)/(rovera+rmaxratio),gamma,rovera)

def getrstarfinaltim(mlossrmax,gamma,rovera=0.1):

    rmaxratio=getrmaxfinal(mlossrmax,gamma,minrmax=True)
    return getrhfinal((rovera+1)/(rovera+rmaxratio),gamma,rovera,minrmax=False)
