import numpy as np

		
x0 = 0.0
y0 = 0.0
z0 = 0.0
x1 = 0.0
y1 = 0.0
z1 = 1.0
N = 6
theta = 180-109
theta = theta*np.pi/180
phi = 120
phi = phi*np.pi/180
Rn = np.matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
 
f = open('coord_6_2.xyz','w')
f.write('{0}\n'.format(N))
f.write('Atoms\n')
f.write('1 {0}\t{1}\t{2}\n'.format(x0, y0, z0))
f.write('1 {0}\t{1}\t{2}\n'.format(x1, y1, z1))

u = np.matrix([[0],[0],[1]])


for i in range(2,N):
    Rz = np.matrix([[np.cos(theta), -np.sin(theta), 0],[np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    #Rx = np.matrix([[1, 0, 0],[0, np.cos(phi), -np.sin(phi)],[0, np.sin(phi), np.cos(phi)]])
    Ry = np.matrix([[np.cos(theta), 0, np.sin(theta)],[0, 1, 0],[-np.sin(theta), 0, np.cos(theta)]])
    Rn = np.matrix(Ry* Rz * Rn)
    u_n = np.matrix(Rn*u)
    #ractual = np.matrix(np.linalg.inv(Rn)*rnew)
    
    x1 = x1 + u_n[0]
    y1 = y1 + u_n[1]
    z1 = z1 + u_n[2]
    f.write('1 ')
    rnew = np.array([x1,y1,z1])
    #f.write('1 {0}\t{1}\t{2}\n'.format(rnew[0], rnew[1], rnew[2]))
    np.savetxt(f, rnew , fmt = '%f', delimiter = '\t',newline = '\t')
    f.write('\n')
    
f.close()
