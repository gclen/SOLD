import gtk
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from ase.gui.widgets import pack
from gettext import gettext as _
from ase.data import chemical_symbols as symbols
from ase.data import atomic_names as names


class Levels(gtk.Window):

    def __init__(self, gui):
        #Set defaults
        self.log_file = None
        self.fchk_file = None
        self.coords = None
        self.scaling_method = 'global'
        
        #Set up window
        gtk.Window.__init__(self)
        self.set_title('Level Diagrams')
        self.gui = gui
        vbox = gtk.VBox()
        vbox.set_border_width(5)

        #Button to select G09 log file
        a = pack(vbox, gtk.Label())
        self.log_entry_box, b = pack(vbox, [gtk.Entry(max=25), 
                                            gtk.Button(_('Choose Gaussian output file'))])
        self.log_entry_box.set_max_length(0)
        self.log_entry_box.connect('activate', self.log_entry)
        b.connect('clicked', self.choose_log_file)

        #Button to select fchk file
        a = pack(vbox, gtk.Label())
        self.fchk_entry_box, b = pack(vbox, [gtk.Entry(max=25), 
                                             gtk.Button(_('Choose fchk file'))])
        self.fchk_entry_box.set_max_length(0)
        self.fchk_entry_box.connect('activate', self.log_entry)
        b.connect('clicked', self.choose_fchk_file)

        #Dial to set number of occupied orbitals
        self.occ_scale = gtk.Adjustment(value=2, lower=0, upper=100, step_incr=1)
        self.occ_spinner = gtk.SpinButton(self.occ_scale, climb_rate=0, digits=0)
        self.occ_spinner.set_update_policy(gtk.UPDATE_IF_VALID)
        self.occ_spinner.set_numeric(True)
        pack(vbox, [gtk.Label(_('Plot ')), self.occ_spinner, gtk.Label(_('occupied states '))])
        self.occ_scale.connect('value-changed', self.scale_occ_orb) 

        #Dial to set number of virtual orbitals
        self.virt_scale = gtk.Adjustment(value=2, lower=0, upper=100, step_incr=1)
        self.virt_spinner = gtk.SpinButton(self.virt_scale, climb_rate=0, digits=0)
        self.virt_spinner.set_update_policy(gtk.UPDATE_IF_VALID)
        self.virt_spinner.set_numeric(True)
        pack(vbox, [gtk.Label(_('Plot ')),
                    self.virt_spinner, gtk.Label(_('unoccupied states '))])
        self.virt_scale.connect('value-changed', self.scale_virt_orb) 

        #Dial to set number of bins
        self.bin_scale = gtk.Adjustment(value=0.3, lower=0.05, upper=0.5, step_incr=0.05)
        self.bin_spinner = gtk.SpinButton(self.bin_scale, climb_rate=0.05, digits=2)
        self.bin_spinner.set_update_policy(gtk.UPDATE_IF_VALID)
        self.bin_spinner.set_numeric(True)
        pack(vbox, [gtk.Label(_(u'Grid spacing (\u212B) ')),
                    self.bin_spinner])
        self.bin_scale.connect('value-changed', self.scale_bin_num) 

        #Button to select atoms
        a = pack(vbox, gtk.Label())
        a, b = pack(vbox, [gtk.Label(_('Select two points and click ')), 
                           gtk.Button(_('Confirm'))])
        b.connect('clicked', self.confirm_points)

        #Button to set scaling
        a = pack(vbox, gtk.Label('Intensity scaling'))
        button = pack(vbox, gtk.RadioButton(None, 'Global'))
        button.connect('toggled', self.set_scaling_method, 'global')
        button.show()
        
        button = pack(vbox, gtk.RadioButton(button, 'Orbital'))
        button.connect('toggled', self.set_scaling_method, 'orbital')
        button.show()

        #Button to set all of the parameters
        a = pack(vbox, gtk.Label())
        a = pack(vbox, gtk.Button(_('Run')))
        a.connect('clicked', self.run_sold)

        # Add elements and show frame
        self.add(vbox)
        vbox.show()
        self.show()

    def confirm_points(self, button):
        indices = np.arange(self.gui.images.natoms)[self.gui.images.selected]
        ordered_indices = self.gui.images.selected_ordered
        n = len(indices)
        self.nselected = n

        if n == 2: 
            self.coords = self.gui.images.P[0][indices]
            self.selected_indices = indices
        else:
            points_error = gtk.MessageDialog(type=gtk.MESSAGE_ERROR, buttons=gtk.BUTTONS_NONE) 
            points_error.set_markup('Please select two points')
            points_error.run()

    def log_entry(self, widget):
        self.log_file = widget.get_text()   

    def fchk_entry(self, widget):
        self.fchk_file = widget.get_text()             

    #Button for getting fchk file
    def choose_fchk_file(self, button):
        chooser = gtk.FileChooserDialog(title=None, action=gtk.FILE_CHOOSER_ACTION_OPEN, 
                       buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        chooser.set_current_folder(os.getcwd())
        
        #Set file filters
        file_filter = gtk.FileFilter()
        file_filter.set_name('Formatted checkpoint files')
        file_filter.add_pattern('*.fchk')
        file_filter.add_pattern('*.fch')
        chooser.add_filter(file_filter)        
        file_filter = gtk.FileFilter()
        file_filter.set_name('All files')
        file_filter.add_pattern('*')
        chooser.add_filter(file_filter)

        response = chooser.run()
        if response == gtk.RESPONSE_OK:
            self.fchk_file = chooser.get_filename()
            self.fchk_entry_box.set_text(self.fchk_file)
        elif response == gtk.RESPONSE_CANCEL:
            print('Closed, no files selected')
        chooser.destroy()     

    #Button for getting gaussian output file
    def choose_log_file(self, button):
        chooser = gtk.FileChooserDialog(title=None, action=gtk.FILE_CHOOSER_ACTION_OPEN, 
                       buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        chooser.set_current_folder(os.getcwd())
        
        #Set file filters
        file_filter = gtk.FileFilter()
        file_filter.set_name('Gaussian output files')
        file_filter.add_pattern('*.log')
        file_filter.add_pattern('*.out')
        chooser.add_filter(file_filter)        
        file_filter = gtk.FileFilter()
        file_filter.set_name('All files')
        file_filter.add_pattern('*')
        chooser.add_filter(file_filter)

        response = chooser.run()
        if response == gtk.RESPONSE_OK:
            self.log_file = chooser.get_filename()
            self.log_entry_box.set_text(self.log_file)
        elif response == gtk.RESPONSE_CANCEL:
            print('Closed, no files selected')
        chooser.destroy()     

    #Button to get desired number of occupied orbitals
    def scale_occ_orb(self, adjustment):
        return True

    #Button to get desired number of virtual orbitals
    def scale_virt_orb(self, adjustment):
        return True

    #Button to get the number of bins    
    def scale_bin_num(self, adjustment):
        return True
    
    #Set the scaling method
    def set_scaling_method(self, widget, data='global'):
        if data == 'orbital' and widget.get_active() == True:
            self.scaling_method = 'orbital'
        elif data == 'global' and widget.get_active() == True:
            self.scaling_method = 'global'

    #Write out the parameters to a config file        
    def write_config(self):     
        occ_num = str(self.occ_spinner.get_value_as_int())
        virt_num = str(self.virt_spinner.get_value_as_int())
        bin_num = str(self.bin_spinner.get_value())
    
        config_file = open('levels_config.txt', 'w')
            
        config_file.write(self.log_file + '\n')
        config_file.write(self.fchk_file + '\n')
        config_file.write(occ_num + '\n')
        config_file.write(virt_num + '\n')
        
        for i in range(len(self.coords)):
            coord_string = str(self.coords[i][0]) + ' ' + str(self.coords[i][1]) + ' ' + str(self.coords[i][2]) 
            config_file.write(coord_string + '\n')

        config_file.write(bin_num + '\n')
        config_file.write(self.scaling_method + '\n')
        config_file.close()

    def plot_levels(self):
        dist = []
        dist_temp = []
        energies = []
        energies_temp = []
        col = []
        col_temp = []

        data_file = open('projected_sums.txt', 'r')

        for line in data_file:
            line_split = line.split()
            
            if line_split[0] == '---':
                dist.append(dist_temp)    
                dist_temp = [] 
                energies.append(energies_temp)
                energies_temp = []
                col.append(col_temp) 
                col_temp = []

            else:
                energies_temp.append(float(line_split[0]))
                dist_temp.append(float(line_split[1]))
                col_temp.append(1-float(line_split[2]))
                  

        data_file.close()

        #Read atoms
        atom_file=open('atom_projections.txt','r')

        atom_points=[]

        for atom in atom_file:
            atom_split=atom.split()
            atom_points.append(float(atom_split[1]))
        
        atom_min = min(atom_points)
        atom_max = max(atom_points)

        dist_adjust = (atom_max+atom_min)/2.0

        dist = [[j-dist_adjust for j in i] for i in dist]

        atom_file.close()      

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        #Adjust x labels
        #Find nearest multiple of 5 lower than the min
        xticks_min = int((atom_min - dist_adjust)- ((atom_min-dist_adjust)%5))
        #Find nearest multiple of 5 higher than the max + 1 to use in range()
        xticks_max = int((atom_max - dist_adjust)+ ((atom_max-dist_adjust)%5)) + 6  

        for i, j in enumerate(dist[0]):
            if j>xticks_max:
                max_index = i    

        for i in range(len(dist)):
            dist[i] = dist[i][:max_index]
            energies[i] = energies[i][:max_index]
            col[i] = col[i][:max_index]

        cdict = {'red': ((0., 1, 1),
                         (0.05, 1, 1),
                         (0.11, 0, 0),
                         (0.66, 1, 1),
                         (0.89, 1, 1),
                         (1, 0.5, 0.5)),
                 'green': ((0., 1, 1),
                           (0.05, 1, 1),
                           (0.11, 0, 0),
                           (0.375, 1, 1),
                           (0.64, 1, 1),
                           (0.91, 0, 0),
                           (1, 0, 0)),
                 'blue': ((0., 1, 1),
                          (0.05, 1, 1),
                          (0.11, 1, 1),
                          (0.34, 1, 1),
                          (0.65, 0, 0),
                          (1, 0, 0))}

        new_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)

        for i in range(len(dist)):
                
            if (i+1 < len(dist)) and abs(energies[i][0]-energies[i+1][0])<0.05:    
                #Get the min and max values of the x axis 
                xmin,xmax=plt.xlim()
                #Add a label indicating the energy shift for degeneracy 
                plt.text(xticks_min-4, energies[i][0],'Shifted for \ndegeneracy')

                for j in range(len(energies[i])):
                    energies[i][j] -= 0.1
                    energies[i+1][j] += 0.1
            

        plt.scatter(dist, energies,cmap=new_cmap, c=col, marker = 's', edgecolor='none')
                            

        plt.xlabel('Distance along v [\AA]')
        plt.ylabel(r'Energy [eV]')          

        #Set the xtics range
        xticks_list=range(xticks_min,xticks_max,5)
        xticks_labels=[str(i) for i in range(xticks_min,xticks_max,5)]
        xticks_labels.append('')

        plt.xticks(xticks_list,xticks_labels)
        
        #Put vertical lines at positions of dotted points
        #Get yrange
        ymin,ymax=plt.ylim()

        #Get desired atom indices
        desired_atom_index = self.selected_indices

        #Store positions in a list
        desired_atoms=[]

        for index in desired_atom_index:
            desired_atoms.append(atom_points[index]-dist_adjust)


        plt.vlines(desired_atoms,ymin,ymax,linestyles='dotted') 

        #Change the xrange
        plt.xlim(((xticks_min-5),(xticks_max+5)))
        #Change the yrange 
        plt.ylim((ymin,ymax))

        plt.show()

    #Button to run the code
    def run_sold(self, button):
        #Check that all parameters have been set
        if self.log_file == None:
            log_error = gtk.MessageDialog(type=gtk.MESSAGE_ERROR, buttons=gtk.BUTTONS_NONE) 
            log_error.set_markup('Please select a Gaussian output file and run again')
            log_error.run()
        if self.fchk_file == None:
            fchk_error = gtk.MessageDialog(type=gtk.MESSAGE_ERROR, buttons=gtk.BUTTONS_NONE) 
            fchk_error.set_markup('Please select a formatted checkpoint file and run again')
            fchk_error.run()
        if self.coords == None:
            points_error = gtk.MessageDialog(type=gtk.MESSAGE_ERROR, buttons=gtk.BUTTONS_NONE) 
            points_error.set_markup('Please select two atoms to use as points')
            points_error.run()            
        
        #Write config and run code
        else:
            self.write_config()
            os.system('SOLD < levels_config.txt')     
            self.plot_levels()
