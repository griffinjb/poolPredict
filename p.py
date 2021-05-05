import numpy as np 
import matplotlib.pyplot as plt 
from numpy.linalg import norm
import pygame as pg 

########################
# Author: 	Josh Griffin
# Date: 	5/3/2021
########################


# distance units in feet

# GLOBALS
BALL_RADIUS = 2.25/12
ROLLING_MU	= .1 		# .005-.015
G 			= 32.174 	# ft/s**2
dt 			= .01 		# second

#sorry
imTBL 	= pg.image.load('table.png')
im8B 	= pg.transform.scale(pg.image.load('8b.png'),[59,59])
imCB 	= pg.transform.scale(pg.image.load('cueball.png'),[59,59])
imSTB 	= pg.transform.scale(pg.image.load('stripe.png'),[59,59])
imSB 	= pg.transform.scale(pg.image.load('solid.png'),[59,59])


### physics equations:

# momentum and therefore velocity is conserved in collsions
# v1**2 = v2a**2 + v2b**2
# post collision v's are orthogonal (if stationary)
# direct hit is full transfer
# force transfer is always c*(p_i - p_j) direction


### sources:
#	https://www.worldscientific.com/doi/pdf/10.1142/9789813276475_0001
#	https://www.mathworks.com/help/stateflow/ug/modeling-the-opening-shot-in-pool.html
#	https://www.real-world-physics-problems.com/physics-of-billiards.html
#	https://www.real-world-physics-problems.com/support-files/cue_ball_traj."%>wde# v1 = v2a + v2b
# 	https://billiards.colostate.edu/faq/physics/physical-properties/
# 	https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-003-signals-and-systems-fall-2011/lecture-videos/MIT6_003F11_lec07.pdf
#  	https://www.researchgate.net/publication/245041330_Collision_of_three_balls_on_a_plane
# 	https://math.stackexchange.com/questions/658871/perfectly-centered-break-of-a-perfectly-aligned-pool-ball-rack/659318#659318
# 	https://github.com/max-kov/pool # Close to what we want
# 	https://github.com/markus-ebke/python-billiards

### basic plan:
# main loop

# while

# 	poll input

	# is clicked
	# mouse position change

# 	compute next state

# 	if time render
###

# Main (temp)
def main():

	N_B = 16
	# init array of balls B = [N_B,4=[px,py,vx,vy]]
	B = np.random.uniform(-1,1,[N_B,4])

# projects u onto v
def proj(u,v):
	return((u@v)/(v@v)*v)

# define get B(+dt)   
# v is user velocity added with strike
def next_B(B,v=np.array([0,0])):
	# get ball interactions

	N_B = B.shape[0]

	B[0,2:] += v

	# Check ball interaction
	Bf = resolve_collision(B)

	# check wall interactions
	# let boundaries be 0,0,4.5,9
	for i in range(N_B):
		if B[i,0] <= 0+BALL_RADIUS:
			Bf[i,0] = 0+BALL_RADIUS+.01
			Bf[i,2] = -B[i,2]
		elif B[i,1] < 0+BALL_RADIUS:
			Bf[i,1] = 0+BALL_RADIUS+.01
			Bf[i,3] = -B[i,3]
		elif B[i,0] > 9-BALL_RADIUS:
			Bf[i,0] = 9-BALL_RADIUS-.01
			Bf[i,2] = -B[i,2]
		elif B[i,1] > 4.5-BALL_RADIUS:
			Bf[i,1] = 4.5-BALL_RADIUS-.01
			Bf[i,3] = -B[i,3]

	# position update and frictional deceleration
	Bf[:,:2] = Bf[:,:2] + Bf[:,2:]*dt
	Bf[:,2:] = Bf[:,2:]*(1- dt*G*ROLLING_MU)

	return(Bf)

def ball_dist(B):
	# upper triangular N_B
	BD = np.zeros(B.shape[0],B.shape[0])

	for i in range(B.shape[0]):
		for j in range(B.shape[0]):
			if j > i:

				BD[i,j] = norm(B[i,:2]-B[j,:2])

# init ball interaction matrix N = [N_B,N_B] ;issym
# define get N from B
def get_interaction(B):
	N = np.zeros([N_B,N_B])
	for i in range(N_B):
		for j in range(i+1,N_B):
			if i != j:
				if np.linalg.norm(B[i,:2]-B[j,:2])<(2*BALL_RADIUS):
					N[i,j] = 1
	return(N)

# 2d
def get_orth_basis(A):
	Af = A.copy()
	Af[:,0] -= proj(A[:,0],A[:,1])
	Af /= norm(Af,axis=0)

def test_triangular_basis():

	c = []
	for i in range(-10,11):
		for j in range(-10,11):
			c.append([i,j])

	c = np.array(c)

	TB = [[1,np.cos(np.pi/3)],[0,np.sin(np.pi/3)]]
	iTB = np.linalg.inv(TB)

	plt.scatter(c[:,0],c[:,1])
	plt.scatter((TB@c.T)[0,:],(TB@c.T)[1,:])
	plt.show()

def collision_2(B):
	# get ball interactions

	Bf = B.copy() 	# Bf - B final - output ball state
	N_B = B.shape[0]

	i = 0
	j = 1

	pi_pj = B[i,:2]-B[j,:2]
	n_pi_pj = norm(pi_pj)
	u_pi_pj = pi_pj/n_pi_pj

	# set post collision velocity
	Bf[i,2:] = B[i,2:] - proj(B[i,2:],u_pi_pj) + proj(B[j,2:],u_pi_pj)
	Bf[j,2:] = B[j,2:] - proj(B[j,2:],u_pi_pj) + proj(B[i,2:],u_pi_pj)

	return(Bf)		

# resolves 3 ball collision
# all balls touching in triangle
def collision_3(B):
	N_B = B.shape[0]
	if N_B != 3:
		print('Warning: More than 3 balls')

	Bf = B.copy()

	P = np.zeros([2,2])

	Bf_pos = [[1,2,0],[0,2,1],[0,1,2]]
	Bf_pos_i = 0

	for i in range(N_B):
		Pi = 0
		for j in range(N_B):

			if i != j:

				# position vectors emanating from i -> j
				P[:,Pi] = B[j,:2] - B[i,:2]
				Pi += 1

		# Inverse triangular basis
		P = P / norm(P,axis=0)

		# basis orthogonalization
		P[:,1] = (P[:,1]+(P[:,1]-P[:,0]))
		P[:,1] /= norm(P[:,1])
		print('corr: ' + str(P[:,0].T@P[:,1]))

		iP = np.linalg.inv(P)

		# invert velocity i 
		vi = iP@B[i,2:].T

		# project 
		vf0 = vi.copy()
		vf1 = vi.copy()
		vf0[1] = 0
		vf1[0] = 0

		vf0 = P@vf0
		vf1 = P@vf1

		# return to triangular basis
			# set velocity final in Bf
		Bf[Bf_pos[Bf_pos_i][0],2:] += vf0
		Bf[Bf_pos[Bf_pos_i][1],2:] += vf1
		Bf[Bf_pos[Bf_pos_i][2],2:] -= (vf0+vf1)
		Bf_pos_i += 1

	return(Bf)

# takes set of intersecting balls
# returns to edge intersection no overlap
# returns backtracked balls (no overlap)
def intersect_backtrack(B):

	# t such that min(max(dist(P*))) == 1+-delta

	# linear approximation
	# P - V*t

	# formulate ball overlap as function of t 

	# P_ = P - V*t_

	# norm(pi-pj - vi*t + vj*t) == 2*BALL_RADIUS

	#  == (2*BALL_RADIUS)**2

	# t = ? 

	# t(vi-vj) = (2*BALL_RADIUS)-DP

	# vi-vj -> vij

	# vij * (2*BALL_RADIUS-||pij||) / ||vij||

	# assuming monotonicity (no jumps)

	# t = (2*BALL_RADIUS-||pij||) / ||vij||

	N_B = B.shape[0]

	DV = np.ones([N_B,N_B])
	DP = 2*BALL_RADIUS*np.ones([N_B,N_B])

	for i in range(N_B):
		for j in range(N_B):
			if j > i:
				DV[i,j] = norm(B[i,2:]-B[j,2:])
				DP[i,j] = norm(B[i,:2]-B[j,:2])

	T = (2*BALL_RADIUS-DP)/DV

	t = np.max(T)

	B[:,:2] -= B[:,2:]*.999999*t

	return(B)

# upper triangular to symmetric matrix
def upp2sym(a):
	return np.where(a,a,a.T)

#										[    ...    ]
# B is array of ball states 			[px,py,vx,vy]
#										[    ...    ]
# returns post collision ball states
def resolve_collision(B):

	# velocity is conserved

	# force is only transferred along contact points

	N_B = B.shape[0]
	B = intersect_backtrack(B)

	# Vi ~ [v,v,v...]
	# Vf ~ [...]
	# ${Vi} == ${Vf}

	# N_B is 2, 3 or 4. 5 is orthogonal -> 0 force transfer -> 2
	
	Bf = B.copy()

	N = get_interaction(B)
	N = upp2sym(N) + np.diag([1.0 for _ in range(N_B)])

	# for col in symmetric interaction matrix
	# 	if 1 -> 2 ball solver
	# 	if 2 -> 3 ball solver
	# 	else undef throw warning, use 2 ball

	for i in range(N_B):
		NC = np.sum(N[i,:])

		if NC == 1:
			_ = None

		elif NC == 2:
			Bc = B[N[i,:]==1,:]
			Bc = collision_2(Bc)
			Bf[N[i,:]==1,:] = Bc

		# elif NC == 3:
		# 	Bc = B[N[i,:]==1,:]
		# 	Bc = collision_3(Bc)
		# 	Bf[N[i,:]==1,:] = Bc
		
		else:
			print('Warning: Poorly resolved collision')
			Bc = B[N[i,:]==1,:]
			Bc = collision_2(Bc)
			Bf[N[i,:]==1,:] = Bc

	return(Bf)

def rack_and_shoot():

	N_B = 16
	B = np.zeros([N_B,4])

	# B[0,:] = (2.0,2.25,30,0)
	B[0,:] = (2.0,2.25,0,0)

	base = np.array([6,2.25])
	loc = base.copy()
	nl = 1

	bi = 1
	for i in range(5):
		loc = base.copy()
		for j in range(nl):
			B[bi,:2] = loc
			bi += 1
			loc += (BALL_RADIUS+.4)*np.array([0,-1])

		nl += 1
		base += (BALL_RADIUS+.4)*np.array([np.cos(np.pi/6),np.sin(np.pi/6)])

	B[:,:2] += np.random.normal(0,.01,[N_B,2])

	return(B)

def plt_balls(B):
	plt.scatter([0,0,9,9],[0,4.5,0,4.5])
	plt.scatter(B[:,0],B[:,1])
	plt.show(block=False)
	plt.pause(.01)
	plt.cla()

def get_game_state():

	# capture screen

	# extract ball position and number of balls

	# position vector cueball - cue
	_=None

def setup_balls():

	# N_B,B,P = get_game_state()

	N_B = 16
	B = rack_and_shoot()

	# init array of balls B = [N_B,4=[px,py,vx,vy]]
	# N_B = 16
	# B = np.random.uniform(0,4.5,[N_B,4])

	# N_B = 4
	# B = np.array([	[1,1,0,0],
	# 				[1,2,0,-10.0],
	# 				[2,1,0,0],
	# 				[2,2,0,-10.0]])

	# N_B = 3
	# B = np.array([
	# 		[2,2.25,10,0],
	# 		[4,2.25+BALL_RADIUS+.0001,0,0],
	# 		[4,2.25-BALL_RADIUS-.0001,0,0]]
	# 	)

	return((B,N_B))

def pxl_to_coord(xy):

	pxpf = (1494 - 85)/9

	cx = (xy[0]-85)/pxpf 

	cy = (xy[1]-84)/pxpf 

	return((cx,cy))	

def coord_to_pxl(xy):
	#85 84 - 1494 787
	# 36**2 
		# -18 -18 px

	# 0,0 -> 85,84

	# px per foot
	pxpf = (1494 - 85)/9
	px = xy[0]*pxpf + 85 - BALL_RADIUS*pxpf

	py = xy[1]*pxpf + 84 - BALL_RADIUS*pxpf

	return((px,py))

def blit_balls(B,gameDisplay):
	# B[0] is cue
	# B[1:-1] is stripe/solid
	# B[-1] is 8ball

	# rescale balls to proper size.

	N_B = B.shape[0]

	for i in range(N_B):

		px,py = coord_to_pxl(B[i,:2])

		if i == 0:
			gameDisplay.blit(imCB, (px,py))
		if i > 0 and i < 8:
			gameDisplay.blit(imSB, (px,py))
		if i > 7 and i < 15:
			gameDisplay.blit(imSTB, (px,py))
		if i == 15:
			gameDisplay.blit(im8B, (px,py))

def render(B,gameDisplay):

	# table -> balls -> guides

	black = (0,0,0)
	white = (255,255,255)

	clock = pg.time.Clock()

	gameDisplay.fill(white)
	gameDisplay.blit(imTBL, (0,0))

	blit_balls(B,gameDisplay)

	pg.display.update()
	clock.tick(60)

def init_display():
	pg.init()
	display_width 	= 1579
	display_height 	= 873
	gameDisplay = pg.display.set_mode((display_width,display_height))
	pg.display.set_caption('8Pool')
	return(gameDisplay)

def UI(B):

	v 	= np.array([0,0])
	gv 	= np.array([0,0])

	FM 	= 30

	mx,my = pxl_to_coord(pg.mouse.get_pos())

	M = np.array([mx,my])

	for event in pg.event.get():

		if event.type == pg.MOUSEBUTTONDOWN:
			if event.button == 1:
				v = B[0,:2] - M

		if event.type == pg.MOUSEBUTTONUP:
			if event.button == 1:
				gv = B[0,:2] - M
		
		if event.type == pg.QUIT:
			pg.quit()

	if np.sum(norm(B[:,2:],axis=1))>1:
		v = np.array([0,0])
	else:
		gv = np.array([0,0])

	return((FM*v,FM*gv))

if __name__ == '__main__':

	B,N_B = setup_balls()

	# B,N_B,P = get_game_state()

	# run sim until collision
	# compute post collision V
	# display V as arrows
	# display transparency cue at collision


	gameDisplay = init_display()

	while 1:

		render(B,gameDisplay)

		v,gv = UI(B)

		B = next_B(B,v)
		B = intersect_backtrack(B)

		# print('v-sum: '+str(np.sum(norm(B[:,2:],axis=1))))



# backtracking
	# may need backtrack 100% 

# 	may need to fix accel backtrack not sure yet.

# vector field regional mask
# perform updates on boundary changes or def 


# # bn ~ a[px,py,vx,vy]
# def F_collision(b0,b1):

# 	# R radius
# 	# kp kv are elasticity constants
# 	# dp = pi - pj center sep balls
# 	# dv = vi-vj is vector diff vel balls
# 	R = BALL_RADIUS
# 	kp,kv = []
# 	dp = b0[:2]-b1[:2]
# 	dv = b0[2:]-b0[2:]
# 	if norm(dp) < 2R:
# 		return(k(2*R-norm(dp))*dp/norm(dp)-k*dv)
# 	else:
# 		return(0)


# class Simulate:
# 	def __init__(self):

# 		self.board = Board()

# 		self.bld = Billiard()

# 		for ball in self.board.balls:

# 			self.bld.add_ball((ball.position[0],ball.position[1]),
# 								(ball.velocity[0],ball.velocity[1]),
# 								radius=BALL_RADIUS)



# 	def advance(self,user_inputs):
# 		# uses ui to determine next state


# class Control:
# 	def __init__(self):
# 		self.mx = 0
# 		self.my = 0
# 		self.isClick = 0

# 	def get(self):
# 		# use pygame to poll UI state


# class Render:
# 	def __init__(self):


# class Guide:
# 	def __init__(self):
# 		# set of guidance lines

# 	def render(self):
# 		# gives bitmap + mask



# class Ball:
# 	def __init__(self,color=1,position=[1,1]):

# 		self.position 		= np.array(position) # [x,y]
# 		self.velocity		= np.array([0,0])
# 		self.color 			= color
# 		self.mass 			= 6 # oz
# 		self.radius 		= 2.25/12 # ft

# 	def render(self):
# 		# gives bitmap + mask



# class Board:
# 	def __init__(self):

# 		self.balls 		= Ball(3)
# 		# self.balls 	= [Ball(0) for _ in range(7)] # solid
# 		# self.balls 	+= [Ball(1) for _ in range(7)] # stripe
# 		# self.balls 	+= [Ball(2)] # eight
# 		# self.balls 	+= [Ball(3)] # cue
# 		self.xmin 	= 0
# 		self.ymin 	= 0
# 		self.xmax 	= 9
# 		self.ymax 	= 4.5

# 		self.setup()

# 	def setup(self):
# 		# arrange balls


# 	def render(self):
# 		# gives bitmap + mask
