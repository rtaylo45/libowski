import glob
import moviepy.editor as mpy

gifName = 'Neutron Precursors'
fps = 10
fileList = glob.glob('*.png')
list.sort(fileList, key=lambda x: float(x.split('_')[1].split('.png')[0]))
clip = mpy.ImageSequenceClip(fileList, fps=fps)
clip.write_gif('{}.gif'.format(gifName), fps=fps)
