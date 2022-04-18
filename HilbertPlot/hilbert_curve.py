import numpy as np

"""convert (x,y) to d"""
def xy2d(n, x, y):
    d = 0
    s = n // 2
    while s > 0:
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * ((3 * rx) ^ ry)
        x, y = rot(n, x, y, rx, ry)
        s //= 2
    return d

"""convert d to (x,y)"""
def d2xy(n, d):
    t = d
    x = 0
    y = 0
    s = 1
    while s < n:
        rx = 1 & (t//2)
        ry = 1 & (t ^ rx)
        x, y = rot(s, x, y, rx, ry)
        x += s * rx
        y += s * ry
        t //= 4
        s *= 2
    return x, y

"""rotate/flip a quadrant appropriately"""
def rot(n, x, y, rx, ry):
    if (ry == 0):
        if (rx == 1):
            x = n-1 - x
            y = n-1 - y

        # Swap x and y
        t = x
        x = y
        y = t
    return x, y

"""3D Hilbert curve"""
def hilbert3(n):
    d1 = np.zeros((n**3), dtype=int)
    d2 = np.zeros((n**3), dtype=int)
    d3 = np.zeros((n**3), dtype=int)
    
    i = 0
    def hilbert_rec(n, x, y, z, dx, dy, dz, dx2, dy2, dz2, dx3, dy3, dz3):
        # http://www.math.uwaterloo.ca/~wgilbert/Research/HilbertCurve/HilbertCurve.html
        # https://stackoverflow.com/a/27080013
        if n == 1:
            nonlocal i, d1, d2, d3
            d1[i] = x
            d2[i] = y
            d3[i] = z
            i += 1
        else:
            n //= 2
            if(dx<0): x-=n*dx
            if(dy<0): y-=n*dy
            if(dz<0): z-=n*dz
            if(dx2<0): x-=n*dx2
            if(dy2<0): y-=n*dy2
            if(dz2<0): z-=n*dz2
            if(dx3<0): x-=n*dx3
            if(dy3<0): y-=n*dy3
            if(dz3<0): z-=n*dz3
            hilbert_rec(n, x, y, z, dx2, dy2, dz2, dx3, dy3, dz3, dx, dy, dz)
            hilbert_rec(n, x+n*dx, y+n*dy, z+n*dz, dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2)
            hilbert_rec(n, x+n*dx+n*dx2, y+n*dy+n*dy2, z+n*dz+n*dz2, dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2)
            hilbert_rec(n, x+n*dx2, y+n*dy2, z+n*dz2, -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3)
            hilbert_rec(n, x+n*dx2+n*dx3, y+n*dy2+n*dy3, z+n*dz2+n*dz3, -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3)
            hilbert_rec(n, x+n*dx+n*dx2+n*dx3, y+n*dy+n*dy2+n*dy3, z+n*dz+n*dz2+n*dz3, -dx3, -dy3, -dz3, dx, dy, dz, -dx2, -dy2, -dz2)
            hilbert_rec(n, x+n*dx+n*dx3, y+n*dy+n*dy3, z+n*dz+n*dz3, -dx3, -dy3, -dz3, dx, dy, dz, -dx2, -dy2, -dz2)
            hilbert_rec(n, x+n*dx3, y+n*dy3, z+n*dz3, dx2, dy2, dz2, -dx3, -dy3, -dz3, -dx, -dy, -dz)
    
    e1 = [1,0,0]
    e2 = [0,1,0]
    e3 = [0,0,1]
    hilbert_rec(n,0,0,0,
             e1[0],e1[1],e1[2],
             e2[0],e2[1],e2[2],
             e3[0],e3[1],e3[2])
    return d1, d2, d3