import controlregion
import antenna
import dlpotential
import boundaryfunctions as bdf
import numpy as np
import regsolve
import sys
import ttk
from Tkinter import *

def compute_phi(cloakframe):
	alpha = cloakframe.alpha.get()
	k = cloakframe.k.get()
	a_rad = cloakframe.a_rad.get()
	n = cloakframe.n.get()
	n_half = cloakframe.n_half.get()
	spacing = cloakframe.spacing.get()
	r1 = cloakframe.r1.get()
	r2 = cloakframe.r2.get()
	theta1 = cloakframe.theta1.get()
	theta2 = cloakframe.theta2.get()
	doffset = cloakframe.doffset.get()
	narcin = cloakframe.narcin.get()
	R = cloakframe.R.get()
	nR = cloakframe.nR.get()
	
	x0 = np.array([1000,0])
	dthetain = (theta2-theta1)/narcin
	narcout = int(np.ceil(r2*(theta2-theta1)/(r1*dthetain)))
	ne = int(np.ceil((r2-r1)/(r1*dthetain)))
	
	a = antenna.AntennaArray(antenna.make_circle_s(a_rad), antenna.circle_ds, spacing, n_half, n)
	c = controlregion.ControlRegion(r1, r2, theta1, theta2, doffset, R, nR, narcin, narcout, ne)

	v1 = bdf.PointSource(k,x0)
	f1 = v1.near_trace(c)
	f2 = np.zeros(c.nR, dtype=complex)
	f = np.concatenate((f1,f2))

	A=dlpotential.compute_k_matrix(a, c, k)
	phi = regsolve.tikhonov_solve2(a, A, f, alpha)
	(Results, Key) = regsolve.analyze(a,c,A,f,phi,k)
	
	L = np.zeros(len(Key), dtype=int)
	for j in range(len(Key)):
		L[j] = len(Key[j])
	s = str(L.max())
	
	output = ""
	format_str = "{:>"+s+"} : {:<6.10f}\n"
	output = output + format_str.format("k", k)
	for entry, datum in zip(Key, np.real(Results)):
		output = output + format_str.format(entry, datum)
	return output
	
class CloakFrame(Frame):
	def __init__(self, parent):
		Frame.__init__(self, parent)
		self.parent = parent
		
		self.alpha = DoubleVar(); self.alpha.set(1e-9);
		self.k = DoubleVar(); self.k.set(1.0);
		self.a_rad = DoubleVar(); self.a_rad.set(0.1);
		self.n = IntVar(); self.n.set(128);
		self.n_half = IntVar(); self.n_half.set(0);
		self.spacing = DoubleVar(); self.spacing.set(1);
		self.r1 = DoubleVar(); self.r1.set(0.11);
		self.r2 = DoubleVar(); self.r2.set(0.14);
		self.theta1 = DoubleVar(); self.theta1.set(-np.pi/4);
		self.theta2 = DoubleVar(); self.theta2.set(np.pi/4);
		self.doffset = DoubleVar(); self.doffset.set(0);
		self.narcin = IntVar(); self.narcin.set(128);
		self.R = DoubleVar(); self.R.set(10);
		self.nR = IntVar(); self.nR.set(256);
		
		self.initUI()
	
	def initUI(self):
		self.grid()
		
		self.alpha_entry = Entry(self, textvariable=self.alpha)
		self.alpha_entry.grid(row=1, column=2, sticky='EW')
		
#		self.alpha_entry.bind("<Return>", self.OnPressEnter)
		self.k_entry = Entry(self, textvariable=self.k)
		self.k_entry.grid(row=2, column=2, sticky='EW')
		self.a_rad_entry = Entry(self, textvariable=self.a_rad)
		self.a_rad_entry.grid(row=3, column=2, sticky='EW')
		self.n_entry = Entry(self, textvariable=self.n)
		self.n_entry.grid(row=4, column=2, sticky='EW')
		self.n_half_entry = Entry(self, textvariable=self.n_half)
		self.n_half_entry.grid(row=5, column=2, sticky='EW')
		self.spacing_entry = Entry(self, textvariable=self.spacing)
		self.spacing_entry.grid(row=6, column=2, sticky='EW')
		self.r1_entry = Entry(self, textvariable=self.r1)
		self.r1_entry.grid(row=7, column=2, sticky='EW')
		self.r2_entry = Entry(self, textvariable=self.r2)
		self.r2_entry.grid(row=8, column=2, sticky='EW')
		self.theta1_entry = Entry(self, textvariable=self.theta1)
		self.theta1_entry.grid(row=9, column=2, sticky='EW')
		self.theta2_entry = Entry(self, textvariable=self.theta2)
		self.theta2_entry.grid(row=10, column=2, sticky='EW')
		self.doffset_entry = Entry(self, textvariable=self.doffset)
		self.doffset_entry.grid(row=11, column=2, sticky='EW')
		self.narcin_entry = Entry(self, textvariable=self.narcin)
		self.narcin_entry.grid(row=12, column=2, sticky='EW')
		self.R_entry = Entry(self, textvariable=self.R)
		self.R_entry.grid(row=13, column=2, sticky='EW')
		self.nR_entry = Entry(self, textvariable=self.nR)
		self.nR_entry.grid(row=14, column=2, sticky='EW')
		
		ttk.Label(self, text="alpha: ").grid(column=1, row=1, sticky=(W,E))
		ttk.Label(self, text="k: ").grid(column=1, row=2, sticky=(W,E))
		ttk.Label(self, text="a: ").grid(column=1, row=3, sticky=(W,E))
		ttk.Label(self, text="n: ").grid(column=1, row=4, sticky=(W,E))
		ttk.Label(self, text="# antenna above/below x-axis: ").grid(column=1, row=5, sticky=(W,E))
		ttk.Label(self, text="antenna spacing: ").grid(column=1, row=6, sticky=(W,E))
		ttk.Label(self, text="control inner radius: ").grid(column=1, row=7, sticky=(W,E))
		ttk.Label(self, text="control outer radius: ").grid(column=1, row=8, sticky=(W,E))
		ttk.Label(self, text="theta_1: ").grid(column=1, row=9, sticky=(W,E))
		ttk.Label(self, text="theta_2: ").grid(column=1, row=10, sticky=(W,E))
		ttk.Label(self, text="d_offset: ").grid(column=1, row=11, sticky=(W,E))
		ttk.Label(self, text="narcin: ").grid(column=1, row=12, sticky=(W,E))
		ttk.Label(self, text="R: ").grid(column=1, row=13, sticky=(W,E))
		ttk.Label(self, text="nR: ").grid(column=1, row=14, sticky=(W,E))
		
		#self.results = StringVar()
		#self.results_label = ttk.Label(self, text=self.results.get()).grid(column=3, row=1, sticky=(W,E))
		#self.results_label['textvariable'] = self.results

		ex_btn = ttk.Button(self,text=u"Compute Phi", command=self.OnCloakButton)
		ex_btn.grid(row=15, column=2)

		self.grid_columnconfigure(0,weight=1)

	def OnCloakButton(self):
		print compute_phi(self)

#	def OnPressEnter(self,event):
#		self.alpha = float(self.alpha_str.get())
#		print "alpha =", self.alpha

def main():
	root = Tk()
	app = CloakFrame(root)
	root.title('Active Field Control Solver')
	root.geometry("400x400")
	root.resizable(True, False)
	root.update_idletasks()
	root.mainloop()
	
if __name__ == "__main__":
	main()
		
#root = Tk()
#root.title("Active Field Control Solver")

#mainframe = ttk.Frame(root, padding="3 3 12 12")
#mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
#mainframe.columnconfigure(0, weight=1)
#mainframe.rowconfigure(0, weight=1)

#alpha = StringVar()
#k = StringVar()
#n = StringVar()
#narcin = StringVar()
#n_antennas = StringVar()
#spacing = StringVar()
#output = StringVar()

#alpha_entry = ttk.Entry(mainframe, width=7, textvariable=alpha)
#alpha_entry.grid(column=0, row=1, sticky=(W,E))
#ttk.Label(mainframe, textvariable=alpha).grid(column=0, row=1, sticky=(W,E))
#ttk.Button(mainframe, text="Compute", command=


