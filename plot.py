
import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.patches import Ellipse, Circle


def plot_emoji(dataset, videoID, tracklist, ax, emotion=None):

	
	ax.plot([1, 9], [5, 5], 'k:', lw=0.5)
	ax.plot([5, 5], [1, 9], 'k:', lw=0.5)

	d = dataset['Emotion'] + '-3' #+ dataset['Level'].map(str)

	val = d.values
	t=tracklist['Title'][videoID-1]
	imgs_name = []

	for i in range(len(val)):
		#if get_emotion_id(val[i]) == emotion:
		imgs_name.append(os.path.join(emoji_path, val[i] + '.png'))
		#imgs_name.append(os.path.join(emoji_path, val[i].split('-')[0] + '-1.png'))

	imgs = []
	for i in range(len(imgs_name)):
		imgs.append(plt.imread(imgs_name[i]))

	dataset['Valence'] = dataset['Valence'].astype(float)
	dataset['Arousal'] = dataset['Arousal'].astype(float)
	x = dataset['Valence'] 
	y = dataset['Arousal'] 

	for x0, y0, path, name in zip(x, y, imgs, imgs_name):
		im = OffsetImage(path, zoom=0.3)
		im.image.axes = ax
		ab = AnnotationBbox(im, (x0, y0), frameon=False)
		ax.add_artist(ab)

	ax.set_xlabel('Valence', fontsize = 16)
	ax.set_ylabel('Arousal', fontsize = 16)
	ax.set_xticks(range(1, 10))
	ax.set_yticks(range(1, 10))
	ax.tick_params(labelsize=14)

	ax.set_xlim(axlim)
	ax.set_ylim(axlim)
	ax.set_aspect('equal', 'box')
	ax.set_title(t, fontsize = 18)
	


###############################     MAIN      #########################################


emoji_path = os.path.abspath(__file__ + '/../emoji/')



fig_dir = os.path.abspath(__file__ + '/../results/')
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)



axlim = [0, 10]

tracklist=pd.read_csv(os.path.abspath(__file__ + '/../videolist.csv'))
tracklist.head()

for i in range(0, 10):

	
	print('Plotting per traccia numero '+str(4*i+1)+'\n')
	print('Plotting per traccia numero '+str(4*i+2)+'\n')
	print('Plotting per traccia numero '+str(4*i+3)+'\n')
	print('Plotting per traccia numero '+str(4*i+4)+'\n')

	dataset1 = pd.read_csv(os.path.abspath(__file__ + '/../input/Results_t'+str(4*i+1)+'.csv'))
	dataset1.head()
	dataset2 = pd.read_csv(os.path.abspath(__file__ + '/../input/Results_t'+str(4*i+2)+'.csv'))
	dataset2.head()
	dataset3 = pd.read_csv(os.path.abspath(__file__ + '/../input/Results_t'+str(4*i+3)+'.csv'))
	dataset3.head()
	dataset4 = pd.read_csv(os.path.abspath(__file__ + '/../input/Results_t'+str(4*i+4)+'.csv'))
	dataset4.head()
	fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(16,9))
	#fig.subplots_adjust(top=0.925, bottom=0.075, wspace=0.2, hspace=0.)
		
	plot_emoji(dataset1, 4*i+1,tracklist, ax[0])
	plot_emoji(dataset2, 4*i+2,tracklist, ax[1])
	plot_emoji(dataset3, 4*i+3,tracklist, ax[2])
	plot_emoji(dataset4, 4*i+4,tracklist, ax[3])

	manager = plt.get_current_fig_manager()
	manager.window.showMaximized()
	plt.tight_layout()

	plt.draw()
	plt.pause(0.1)
	fig.savefig(os.path.join(fig_dir, 'grafici_t'+str(4*i+1)+'-'+str(4*i+4)+'.jpg' ))
	plt.close(fig)


# Clear the __pycache__ folder
pycache_path = os.path.abspath(__file__ + '/../__pycache__')
if os.path.exists(pycache_path):
	import shutil
	shutil.rmtree(pycache_path)