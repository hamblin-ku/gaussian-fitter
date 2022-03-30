import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
from scipy.optimize import curve_fit
import numpy as np
from scipy import asarray as ar,exp
import tkinter as tk
#import tkFileDialog
from tkinter import ttk
from tkinter import *
from tkinter.ttk import *
#from tkinter.ttk import Frame, Button, Style, Label
from astropy.io import ascii
from astropy.table import Table, Column

params = {
    'text.latex.preamble': '\\usepackage{gensymb}',
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray', #gray
    'axes.grid': False,
    'savefig.dpi': 150,  # to adjust notebook inline plot size
    'axes.labelsize': 12, # fontsize for x and y labels (was 10); 40
    'axes.titlesize': 20, #40
    'font.size': 28, # was 10 #28
    'legend.fontsize': 30, # was 10; 28
    'xtick.labelsize': 10, #28
    'ytick.labelsize': 10, #28
    'text.usetex': True,
    'figure.autolayout': True,
    #'figure.figsize': [3.39, 2.10],
    'font.family': 'serif',
}
matplotlib.rcParams.update(params)

f = Figure(figsize=(5,5), dpi=100)
a = f.add_subplot(111)
a.set_xlabel(r'Wavelength\ (nm)')
a.set_ylabel(r'Intensity')

#fig, ax = plt.subplots()
def gaussian(x, height, center, sigma, offset):
    return height*np.exp(-(x - center)**2/(2*sigma**2)) + offset
    
def two_gaussians(x, h1, c1, w1, h2, c2, w2, offset):
    return (gaussian(x, h1, c1, w1, offset=0) +
            gaussian(x, h2, c2, w2, offset=0) + offset)

class GaussFit(Frame):
    def __init__(self):
        Frame.__init__(self)
        self.initUI()

    def initUI(self):
        textpadding =1
        
        ## Params
        self.N_gauss = 1
        self.fileName = ""
        self.data = Table(data = [[0,0],[0,0] ])
        self.notLoaded = True
        ##
        
        self.style = Style()
        self.style.theme_use("aqua") #default
        self.master.title("Guassian Fitter")
        self.pack(fill=BOTH, expand = 1)
        
        ## Frame containers
        mainControls = tk.Frame(self, background = "#3d3c3c")
        mainControls.pack(ipady=10,ipadx=10, fill=BOTH, expand  = True, side = LEFT)
       
        self.canvasPlot = FigureCanvasTkAgg(f, self) #, width = 500, height = 1000)
        self.canvasPlot.get_tk_widget().pack(fill=BOTH, expand = True, side = LEFT)
        
        ###
        
        ## Content for Frames
        
        ###: Frames for # of Gaussians
        selectFrame = tk.Frame(mainControls, background = "#3d3c3c")
        selectFrame.pack(fill=Y, expand = False, side = TOP) #Fill along Y, dont take up whole frame, start at top

        ###
        selectionLbl = Label(selectFrame, text=" Number of Gaussians to fit: ", font=("Helvetica", 12))
        selectionLbl.pack(side=LEFT, expand=True, fill=X)
        
        self.selectionBox = Spinbox(selectFrame, justify = 'center' ,from_=1, to=2, width=3, command=self.set_N)
        self.selectionBox.pack(side = RIGHT, expand=False)
        self.selectionBox.insert(0,1)
        ###
    
        
        
        ###: Data load
        dataButton = Button(mainControls,text = "Load File", command = self.loadFile)
        dataButton.pack(fill=X, pady=0)
        ###
        
        ###: Perform the fit:
        fitButton = Button(mainControls,text = "Perform Fit", command = self.fit)
        fitButton.pack(fill =X, pady=0)
        ###
        
        ###: Refresh
        refreshButton = Button(mainControls, text = "Refresh", command = self.refreshFigure)
        refreshButton.pack(fill = X, pady=0)
        ###
        
        ###: Save the plot:
        SaveButton = Button(mainControls, text = "Save Plot", command = self.savePlot)
        SaveButton.pack(fill =X, pady=0)
        ###
        
        
        ### Coordinate limit Entries
        limitFrameX = tk.Frame(mainControls, background = "#3d5c3c")
        limitFrameX.pack(fill = Y, expand = False, side = TOP, pady = (10,0))
        x_limits = Label(limitFrameX, text = " X limits ")
        x_limits.pack( expand = True, )
        self.lower_xBox = Entry(limitFrameX, width = 9, justify = 'center')#, exportselection = 0)
        self.lower_xBox.pack(side = LEFT, expand = False)
        self.higher_xBox = Entry(limitFrameX, justify = 'center', exportselection = 0,  width = 9)
        self.higher_xBox.pack(side = RIGHT, expand = False)
    
        limitFrameY = tk.Frame(mainControls, background = "#3d5c3c")
        limitFrameY.pack(fill = Y, expand = False, side = TOP)
        y_limits = Label(limitFrameY, text = " Y limits ")
        y_limits.pack( expand = True)
        self.lower_yBox = Entry(limitFrameY, justify = 'center', exportselection = 0,  width = 9)
        self.lower_yBox.pack(side = LEFT, expand = False)
        self.higher_yBox = Entry(limitFrameY, justify = 'center', exportselection = 0,  width = 9)
        self.higher_yBox.pack(side = RIGHT, expand = False)
        ###

        ###
        centerFrame = tk.Frame(mainControls,background = "#3d3c3c")
        centerFrame.pack(fill = Y, expand = False, side = TOP, pady = (10,0))
        GuessParams = Label(centerFrame, text = " Fit Parameter Guesses ", anchor = CENTER,width = 16)
        GuessParams.pack(fill = X,side = TOP)
        self.c1_Label = Label(centerFrame, text = "c_1", width = 3)
        self.c1_Label.pack(fill = Y,side = LEFT)
        self.c1_box = Entry(centerFrame, justify = 'center', exportselection = 0, width = 5)
        self.c1_box.pack(side = LEFT, expand = False)
        self.c2_Label = Label(centerFrame, text = "c_2", width = 3)
        self.c2_Label.pack(fill = Y, side = LEFT)
        self.c2_box = Entry(centerFrame, justify = 'center', exportselection = 0, width = 5)
        self.c2_box.pack(side = LEFT, expand = False)
        ###
                
        ###
        heightFrame = tk.Frame(mainControls,background = "#3d5c3c")
        heightFrame.pack(fill = Y, expand = False, side = TOP)
        self.h1_Label = Label(heightFrame, text = "h_1", width = 3)
        self.h1_Label.pack(fill = Y,side = LEFT)
        self.h1_box = Entry(heightFrame, justify = 'center', exportselection = 0, width = 5)
        self.h1_box.pack(side = LEFT, expand = False)
        self.h2_Label = Label(heightFrame, text = "h_2", width = 3)
        self.h2_Label.pack(fill = Y, side = LEFT)
        self.h2_box = Entry(heightFrame, justify = 'center', exportselection = 0, width = 5)
        self.h2_box.pack(side = LEFT, expand = False)
        
        self.h2_box.configure(state = 'disabled') #disable h2 and c2 since by default N_gauss = 1
        self.c2_box.configure(state = 'disabled')
        ###
        
        ###
        yLabelFrame = tk.Frame(mainControls,background = "#3d3c3c")
        yLabelFrame.pack(fill = Y, expand = False, side = TOP, pady = (10,0))
        yLabelHeader = Label(yLabelFrame, text = "y-axis Label text", anchor = CENTER, width = 20)
        yLabelHeader.pack(fill = X,side = TOP)
        self.ylabel_box = Entry(yLabelFrame, justify = 'center', exportselection = 0, width = 20)
        self.ylabel_box.pack(side = BOTTOM, expand = False)
        self.ylabel_box.insert(0, "Intensity")
        ###
        
        ###
        xLabelFrame = tk.Frame(mainControls,background = "#3d3c3c")
        xLabelFrame.pack(fill = Y, expand = False, side = TOP, pady = (0,0))
        xLabelHeader = Label(xLabelFrame, text = "x-axis Label text", anchor = CENTER, width = 20)
        xLabelHeader.pack(fill = X,side = TOP)
        self.xlabel_box = Entry(xLabelFrame, justify = 'center', exportselection = 0, width = 20)
        self.xlabel_box.pack(side = BOTTOM, expand = False)
        self.xlabel_box.insert(0, "Wavelength (nm)")
        ###
        
    def set_N(self):
        val = int(self.selectionBox.get())
        if self.N_gauss == 1 and val == 2:
            self.h2_box.configure(state = 'normal')
            self.c2_box.configure(state = 'normal')
        elif self.N_gauss == 2 and val == 1:
            self.h2_box.delete(0, END)
            self.c2_box.delete(0, END)
            self.h2_box.configure(state = 'disabled')
            self.c2_box.configure(state = 'disabled')
        self.N_gauss = val
        return

    def clearFit(self):
        if self.N_gauss == 1:
            self.h1_box.delete(0, END)
            self.c1_box.delete(0, END)
        if self.N_gauss == 2:
            self.h1_box.delete(0, END)
            self.h2_box.delete(0, END)
            self.c1_box.delete(0, END)
            self.c2_box.delete(0, END)
    

    def loadFile(self):
        self.fileName = tk.filedialog.askopenfilename(initialdir = "/", title = "Select Data File", filetypes = (("csv files", "*.csv"),("txt files", "*.txt")))
        self.data = Table(ascii.read(self.fileName))
        self.notLoaded = False
        self.clearFit()

        self.plot()
    
    def savePlot(self):
        save_name = tk.filedialog.asksaveasfilename(initialfile = "untitled_spectra.pdf", defaultextension=".pdf", filetypes=[("All Files","*.*"),("Text Documents","*.txt")])
        if save_name is None:
            return
        f.savefig(save_name)
        
    def refreshFigure(self):
        
        a.set_ylabel(r''+self.ylabel_box.get())
        a.set_xlabel(r''+self.xlabel_box.get())
        
        if self.notLoaded == True:
            #print "######## \n \n",self.lower_xBox.get()+1, "############## \n \n"
            a.set_xlim(float(self.lower_xBox.get()), float(self.higher_xBox.get()))
            #print "######## \n \n",self.lower_xBox.get(), "############## \n \n"
            a.set_ylim(float(self.lower_yBox.get()), float(self.higher_yBox.get()))
        elif self.notLoaded == False:
            x = self.data['col1']
            y = self.data['col2']
            x_min, x_max, y_min, y_max = min(x), max(x), min(y), max(y)
            a.set_xlim([x_min, x_max])
            a.set_ylim([y_min, y_max])
        
            self.lower_xBox.delete(0, END)
            self.lower_yBox.delete(0, END)
            self.higher_xBox.delete(0, END)
            self.higher_yBox.delete(0, END)
        
            self.lower_xBox.insert(END, x_min)
            self.lower_yBox.insert(END, y_min)
            self.higher_xBox.insert(END, x_max)
            self.higher_yBox.insert(END, y_max)
            self.notLoaded = True
            
            
            
        self.canvasPlot.draw()
        return
    
    def clearFig(self):
        a.clear()
        self.plot()
    
    def getFileName(self):
        return self.fileName
    
    def plot(self):
        #self.data = Table(ascii.read(self.getFileName()))
        a.clear()
        a.set_ylabel(r''+self.ylabel_box.get())
        a.set_xlabel(r''+self.xlabel_box.get())
        
        a.scatter(self.data['col1'], self.data['col2'], c = 'k', marker = '.') ##Figure out how to embed matplotlib plot into a plotFrame
        index = len(self.fileName) - 1
        while self.fileName[index] != "/":
            index -= 1
        name_pre = self.fileName[(index+1):]
        name = ""
        for c in name_pre:
            if c not in ["_"]:
                name += c
        a.set_title(r''+str(name))
        
        self.refreshFigure()
        return

    def fit(self):
        self.clearFig()
        x = self.data['col1']
        y = self.data['col2']
        
        if self.N_gauss == 1:
            h_1 = float(self.h1_box.get())
            c_1 = float(self.c1_box.get())
            guess = [h_1, c_1, .1, 0]
    
            popt, pcov = curve_fit(gaussian, x, y, p0 = guess)
            err = np.sqrt(np.diag(pcov))

            print('\n********************')
            print('***Fit Parameters***')
            print('********************\n ')
            print( f'height = {popt[0]:.3f} +- {err[0]:.2f}')
            print( f'center = {popt[1]:.3f} +- {err[1]:.2f}')
            print( f'sigma = {popt[2]:.3f} +- {err[2]:.2f}\n')
            print('********************')
            
            x_plot = np.linspace(x[0], x[-1], len(x)*10)
            y_interp = gaussian(x_plot, *popt)
            a.plot(x_plot, y_interp, c = 'r')
            self.refreshFigure()


        elif self.N_gauss == 2:
            h_1 = float(self.h1_box.get())
            h_2 = float(self.h2_box.get())
            c_1 = float(self.c1_box.get())
            c_2 = float(self.c2_box.get())
            
            guess = [h_1, c_1, .1, h_2, c_2,.1, 0]
            
            popt, pcov = curve_fit(two_gaussians, x, y, p0 = guess)
            err = np.sqrt(np.diag(pcov))

            print('\n********************')
            print('***Fit Parameters***')
            print('********************\n ')
            print( f'height_1 = {popt[0]:.3f} +- {err[0]:.2f}')
            print( f'center_1 = {popt[1]:.3f} +- {err[1]:.2f}')
            print( f'sigma_1 = {popt[2]:.3f} +- {err[2]:.2f}\n')

            print( f'height_2 = {popt[3]:.3f} +- {err[3]:.2f}')
            print( f'center_2 = {popt[4]:.3f} +- {err[4]:.2f}')
            print( f'sigma_2 = {popt[5]:.3f} +- {err[5]:.2f}\n')
            print('********************')
            
            x_plot = np.linspace(x[0], x[-1], len(x)*10)
            y_interp = two_gaussians(x_plot, *popt)
            a.plot(x_plot, y_interp, c = 'r')
            self.refreshFigure()



class StartPage(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        label = tk.Label(self, text="Start Page", font=LARGE_FONT)
        label.pack(pady=10,padx=10)
        
        button = ttk.Button(self, text="One Gaussian Peak",command=lambda:controller.show_frame(PageOne))
        button.pack()
        
        button2 = ttk.Button(self, text="Two Gaussian Peak",command=lambda: controller.show_frame(PageTwo))
        button2.pack()
        
        #button3 = ttk.Button(self, text="Graph Page", command=lambda: controller.show_frame(PageThree))
        #button3.pack()
class PageOne(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Single Gaussian", font=LARGE_FONT)
        label.pack(pady=10,padx=10)
        
        button1 = ttk.Button(self, text="Back to Home",command=lambda: controller.show_frame(StartPage))
        button1.pack()
     
        button2 = ttk.Button(self, text="Page Two",command=lambda: controller.show_frame(PageTwo))
        button2.pack()

class PageTwo(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Two Gaussians", font=LARGE_FONT)
        label.pack(pady=10,padx=10)
        
        button1 = ttk.Button(self, text="Back to Home",command=lambda: controller.show_frame(StartPage))
        button1.pack()
        
        button2 = ttk.Button(self, text="Single Gaussian",command=lambda: controller.show_frame(PageOne))
        button2.pack()
         
        canvas = FigureCanvasTkAgg(f, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
         
        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

def main():
    
    root = tk.Tk()
    app = GaussFit()
    root.mainloop()


if __name__ == '__main__':
    main()
