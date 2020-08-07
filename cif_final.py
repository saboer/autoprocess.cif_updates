from .base import Base
from jinja2 import Environment, FileSystemLoader
from processing.models import Collection, Processing, setup
from beamline import variables as blconfig
import re 
import pandas as pd
import numpy as np
import math 


def run_test():
    beamline = 'MX2' # beamline in collection objecxt
    energy_in_kev = 13.25
    detector = 'eiger' #detector_type in collection object
    cryojet_temperature = 110.2
    project_dir = '.'
    write_template_file(project_dir, beamline, detector, energy_in_kev, cryojet_temperature)

def write_template_file(project_dir, beamline, detector, energy_in_kev, cryojet_temperature, crystal_in_monochromator, sample_desc, sample_colour, sample_size_min, sample_size_mid, sample_size_max):
    if beamline == 'MX1':
        beamline_text = 'MX1 Beamline Australian Synchrotron'
    elif beamline == 'MX2':
        beamline_text = 'MX2 Beamline Australian Synchrotron'
    else:
        raise Exception('Unknown beamline')

    if beamline == 'MX1':
        template_mx = '%s/%s' % (project_dir, '/autoprocess_new.cif_mx1')
        template_mx_sadw = '%s/%s' % (project_dir, 'sadabs_w/autoprocess_sadabs_w.cif_mx1')
        template_mx_sadm = '%s/%s' % (project_dir, 'sadabs_m/autoprocess_sadabs_m.cif_mx1')
        template_mx_sads = '%s/%s' % (project_dir, 'sadabs_s/autoprocess_sadabs_s.cif_mx1')
    elif beamline == 'MX2':
        template_mx = '%s/%s' % (project_dir, '/autoprocess_new.cif_mx2')
        template_mx_sadw = '%s/%s' % (project_dir, 'sadabs_w/autoprocess_sadabs_w.cif_mx2')
        template_mx_sadm = '%s/%s' % (project_dir, 'sadabs_m/autoprocess_sadabs_m.cif_mx2')
        template_mx_sads = '%s/%s' % (project_dir, 'sadabs_s/autoprocess_sadabs_s.cif_mx2')
    else:
        raise Exception('Unknown beamline')

    if detector == 'eiger' and beamline == 'MX1':
        detector_text = 'Dectris Eiger2 9M'
    elif detector == 'eiger' and beamline == 'MX2':
        detector_text = 'Dectris Eiger 16M'
    else:
        raise Exception('Unknown detector')

    env = Environment(
        loader=FileSystemLoader('/xray/software/Python/libraries/mx_auto_dataset/templates'))

    template = env.get_template('cx_template.cif')
    KEV_TO_ANGSTROM = 12.398420
    wavelength = KEV_TO_ANGSTROM/float(energy_in_kev)

    if crystal_in_monochromator == 'DC':
        crystal_in_monochromator = 'Silicon Double Crystal'
    elif crystal_in_monochromator == 'CC':
        crystal_in_monochromator = 'Silicon Channel Cut Crystal'
    else:
        raise Exception('Unknown monochromator type')
    with open('%s/%s' % (project_dir, 'IDXREF.LP'), 'r') as index_file:
        contents = index_file.read()
    x = re.search('AUTOINDEXING IS BASED ON',contents)
    y = re.search('OUT OF',contents)    
    index_refs = contents[x.end():y.start()].strip(' ')

    theta_min, theta_max = index_angles(project_dir)
    theta_min = f'{theta_min:.2f}'
    theta_max = f'{theta_max:.2f}'
    temperature = f'{cryojet_temperature:.0f}' 
    size_min = sample_size_min / 1000
    size_mid = sample_size_mid / 1000
    size_max = sample_size_max / 1000
    
    Tmin_w, Tmax_w = abs_Tminmax(sad_w)
    Tmin_m, Tmax_m = abs_Tminmax(sad_m)
    Tmin_s, Tmax_s = abs_Tminmax(sad_s)    

    with open(template_mx, 'w') as template_file:
        template_file.write(template.render(detector=detector_text, beamline=beamline_text, 
                                            wavelength='%.6f' % wavelength, index = index_refs, 
                                            temperature=temperature, crystal=crystal_in_monochromator, 
                                            description = sample_desc, colour = sample_colour, 
                                            size_min = size_min, size_mid = size_mid, 
                                            size_max = size_max, theta_min = theta_min, theta_max = theta_max, 
                                            Tmin = 'Value not reported by XDS', Tmax = 'Value not reported by XDS', 
                                            abs = 'XDS (Kabsch, 2010)'))

    with open(template_mx_sadw, 'w') as template_file:
        template_file.write(template.render(detector=detector_text, beamline=beamline_text, 
                                            wavelength='%.6f' % wavelength, index = index_refs, 
                                            temperature=temperature, crystal=crystal_in_monochromator, 
                                            description = sample_desc, colour = sample_colour, 
                                            size_min = size_min, size_mid = size_mid, 
                                            size_max = size_max, theta_min = theta_min, theta_max = theta_max, 
                                            Tmin = Tmin_w, Tmax = Tmax_w, abs = 'sadabs (Bruker, 2001)'))

    with open(template_mx_sadm, 'w') as template_file:
        template_file.write(template.render(detector=detector_text, beamline=beamline_text, 
                                            wavelength='%.6f' % wavelength, index = index_refs, 
                                            temperature=temperature, crystal=crystal_in_monochromator, 
                                            description = sample_desc, colour = sample_colour, 
                                            size_min = size_min, size_mid = size_mid, 
                                            size_max = size_max, theta_min = theta_min, theta_max = theta_max, 
                                            Tmin = Tmin_m, Tmax = Tmax_m, abs = 'sadabs (Bruker, 2001)'))
    
    with open(template_mx_sads, 'w') as template_file:
        template_file.write(template.render(detector=detector_text, beamline=beamline_text, 
                                            wavelength='%.6f' % wavelength, index = index_refs, 
                                            temperature=temperature, crystal=crystal_in_monochromator, 
                                            description = sample_desc, colour = sample_colour, 
                                            size_min = size_min, size_mid = size_mid, 
                                            size_max = size_max, theta_min = theta_min, theta_max = theta_max, 
                                            Tmin = Tmin_s, Tmax = Tmax_s, abs = 'sadabs (Bruker, 2001)'))
        
def index_angles(project_dir):
    with open('%s/%s' % (project_dir, 'IDXREF.LP')) as f:
        contents = f.read()
    w = re.search('DETECTOR_DISTANCE=', contents)
    x = re.search('ORGX=', contents)
    y = re.search('ORGY=', contents)
    z = re.search('NUMBER', contents)
    det_mm = float(contents[w.end():x.start()].strip(' '))
    ORGX = float(contents[x.end():y.start()].strip(' '))
    ORGY = float(contents[y.end():z.start()].strip(' '))

    df = pd.read_csv('%s/%s' % (project_dir, 'SPOT.XDS'), sep='\s+', names=list('XYPchkl')) 
    df = df[(df[["h", "k", "l"]].T != 0).any()]
    x = df['X'] - ORGX
    y = df['Y'] - ORGY
    theta_max = math.degrees(math.atan((max(np.hypot(x, y))*0.075) / det_mm))/2
    theta_min = math.degrees(math.atan((min(np.hypot(x, y))*0.075) / det_mm))/2
    
    return (theta_min, theta_max)

def abs_Tminmax(sad):
    with open('%s/%s' % (sad, 'sad.abs')) as g:
        contents = g.read()
    a = re.search('maximum transmission:', contents)
    b = re.search('The ratio', contents)
    Tminmax = contents[a.end():b.start()].split(' ')
    Tmin = float(Tminmax[2])
    Tmax = float(str.rstrip(Tminmax[4]))
    return (Tmin, Tmax)

class Cif(Base):

    def __init__(self, run_name, *args, **kwargs):
        super(Cif, self).__init__()

    def process(self, **kwargs):
        if kwargs['collection_id']:
            coll = Collection(kwargs['collection_id'])
        else:
            setup(blconfig.get_database())
            proc = Processing(kwargs['dataset_id'])
            coll = Collection(str(proc.collection_id.id))

        try:
            cryo_temp = coll.cryo_temperature
        except AttributeError:
            cryo_temp = None
        try:
            crystal_in_monochromator = coll.crystal_in_monochromator
        except AttributeError:
            crystal_in_monochromator = None
        try:
            sample_desc = coll.sample_desc
        except AttributeError:
            sample_desc = 'user input'
        try:
            sample_colour = coll.sample_colour
        except AttributeError:
            sample_colour = 'user input'
        try:
            sample_size_min = coll.sample_size_min
        except AttributeError:
            sample_size_min = 'user input'
        try:
            sample_size_mid = coll.sample_size_mid
        except AttributeError:
             sample_size_mid = 'user input'
        try:
            sample_size_max = coll.sample_size_max
        except AttributeError:
             sample_size_max = 'user input'

        write_template_file(self.project_dir, coll.beamline, coll.detector_type, coll.energy, cryo_temp, crystal_in_monochromator, sample_desc, sample_colour, sample_size_min, sample_size_mid, sample_size_max)

if __name__ == '__main__':
    run_test()
