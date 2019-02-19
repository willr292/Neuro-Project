from neuron import h, gui
from matplotlib import pyplot
import numpy

soma = h.Section(name='soma')
dend = h.Section(name='dend')

h.psection(sec=soma)

dend.connect(soma(1))

h.psection(sec=dend)

# Surface area of cylinder is 2*pi*r*h (sealed ends are implicit).
# Here we make a square cylinder in that the diameter
# is equal to the height, so diam = h. ==> Area = 4*pi*r^2
# We want a soma of 500 microns squared:
# r^2 = 500/(4*pi) ==> r = 6.2078, diam = 12.6157
soma.L = soma.diam = 12.6157 # Makes a soma of 500 microns squared.
dend.L = 200 # microns
dend.diam = 1 # microns
print("Surface area of soma = {}".format(soma(0.5).area()))
shape_window = h.PlotShape()
shape_window.exec_menu('Show Diam')

for sec in h.allsec():
    sec.Ra = 100    # Axial resistance in Ohm * cm
    sec.cm = 1      # Membrane capacitance in micro Farads / cm^2

# Insert active Hodgkin-Huxley current in the soma
soma.insert('hh')
for seg in soma:
    seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
    seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
    seg.hh.gl = 0.0003    # Leak conductance in S/cm2
    seg.hh.el = -54.3     # Reversal potential in mV

# Insert passive current in the dendrite
dend.insert('pas')
for seg in dend:
    seg.pas.g = 0.001  # Passive conductance in S/cm2
    seg.pas.e = -65    # Leak reversal potential mV

stim = h.IClamp(dend(1))

stim.delay = 5
stim.dur = 1
stim.amp = 0.1

v_vec = h.Vector()        # Membrane potential vector
t_vec = h.Vector()        # Time stamp vector
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)
simdur = 25.0

h.tstop = simdur
h.run()

pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
pyplot.plot(t_vec, v_vec)
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()

pyplot.figure(figsize=(8,4))
step = 0.075
num_steps = 4
for i in numpy.linspace(step, step*num_steps, num_steps):
    stim.amp = i
    h.tstop = simdur
    h.run()
    pyplot.plot(t_vec, v_vec, color='black')

pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()

dend_v_vec = h.Vector()        # Membrane potential vector
dend_v_vec.record(dend(0.5)._ref_v)

pyplot.figure(figsize=(8,4))
for i in numpy.linspace(step, step*num_steps, num_steps):
    stim.amp = i
    h.tstop = simdur
    h.run()
    soma_plot = pyplot.plot(t_vec, v_vec, color='black')
    dend_plot = pyplot.plot(t_vec, dend_v_vec, color='red')

# After looping, actually draw the image with show.
# For legend labels, use the last instances we plotted
pyplot.legend(soma_plot + dend_plot, ['soma', 'dend'])
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()