import time
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import ZeroPadding2D, GaussianDropout, Dense, Dropout, Activation, Reshape, Flatten, GaussianNoise
from tensorflow.keras.layers import Lambda, Add, Input,Concatenate, Cropping2D, UpSampling2D, ReLU, SpatialDropout2D, LeakyReLU,Conv3D, SeparableConv2D, Conv2D, Conv2DTranspose, BatchNormalization, AveragePooling2D, MaxPooling2D
from tensorflow.keras.optimizers import RMSprop, Adam
from tensorflow.keras import regularizers
from tensorflow.keras.callbacks import ModelCheckpoint
import h5py
import numpy as np
from copy import deepcopy
from scipy.ndimage import gaussian_filter







def pixel_shuffle(scale):
	return lambda x: tf.nn.depth_to_space(x,scale)

def res_block(x_in,filters,scaling):
	x = Conv2D(filters,3,padding='same',activation='relu')(x_in)
	x = Conv2D(filters,3,padding='same',activation='relu')(x)
	x = Lambda(lambda t: t*0.1)(x)
	x = Add()([x_in,x])
	return x

def upsample(x,scale,num_filters):
	def upsample_1(x,factor,**kwargs):
		x = Conv2D(num_filters*(factor**2),3,padding='same',**kwargs)(x)
		return Lambda(pixel_shuffle(scale=factor))(x)
	if scale == 2:
		x = upsample_1(x,2,name='conv2d_1_scale_2')
	return x

echo_in = Input(shape=(256,256,2))
l1 = Conv2D(activation=None,filters=64,kernel_size=(3,3),padding='same',data_format='channels_last',name='conv2d_vgg_0',trainable=True)(echo_in)
l1a = LeakyReLU()(l1)
l2 = MaxPooling2D(pool_size=(2,2))(l1a)

l3 = Conv2D(filters=128,kernel_size=(3,3),activation=None,padding='same',name='conv2d_vgg_1',trainable=True)(l2)
l3a = LeakyReLU()(l3)
l4 = MaxPooling2D(pool_size=(2,2))(l3a)

l5 = Conv2D(filters=256,kernel_size=(3,3),activation=None,padding='same',name='conv2d_vgg_2',trainable=True)(l4)
l5a = LeakyReLU()(l5)
l6 = Conv2D(filters=256,kernel_size=(3,3),activation=None,padding='same',name='conv2d_vgg_3',trainable=True)(l5a)
l6a = LeakyReLU()(l6)
l7 = MaxPooling2D(pool_size=(2,2))(l6a)

l8 = Conv2D(filters=512,kernel_size=(3,3),activation=None,padding='same',name='conv2d_vgg_4',trainable=True)(l7)
l8a = LeakyReLU()(l8)
l9 = Conv2D(filters=512,kernel_size=(3,3),activation=None,padding='same',name='conv2d_vgg_5',trainable=True)(l8a)
l9a = LeakyReLU()(l9)

l10 = MaxPooling2D(pool_size=(2,2))(l9a)
l11 = Conv2D(filters=512,kernel_size=(3,3),activation=None,padding='same',name='conv2d_vgg_6',trainable=True)(l10)
l11a = LeakyReLU()(l11)
l12 = Conv2D(filters=512,kernel_size=(3,3),activation=None,padding='same',name='conv2d_vgg_7',trainable=True)(l11a)
l12a = LeakyReLU()(l12)
l13 = MaxPooling2D(pool_size=(2,2))(l12a)
l14 = SpatialDropout2D(0.5)(l13)
l15 = Flatten()(l14)
l16 = Dense(19*106,activation=None,name='dense_vgg_1',trainable=True)(l15)
l16a = LeakyReLU()(l16)
l17 = Dense(19*106,activation=None,name='dense_vgg_2',trainable=True)(l16a)
l18 = LeakyReLU()(l17)

x_in = Reshape((19,106,1))(l18)

x = b = Conv2D(64,3,padding='same')(x_in)

for i in range(16):
	b = res_block(b,64,0.1)
b = Conv2D(64,3,padding='same')(b)
x = Add()([x,b])

x = upsample(x,2,64)

xout = Conv2D(1,3,padding='same')(x)

rms = Adam(lr=1e-4)
mymodel = Model(inputs=[echo_in],outputs=[xout])
mymodel.compile(loss=['mean_absolute_error'],optimizer=rms,metrics=['mean_squared_error'])

print('ready to goo!!!')
mymodel.summary()

data_tr = np.load('/home/echo_hr_108_train.npy')
data_val = np.load('/home/echo_hr_108_val.npy')
data_mean = np.load('/home/echo_hr_108_mean.npy')
data_tr = (data_tr-data_mean)
data_val = (data_val-data_mean)
ims = np.load('/home/sar_1015_lr.npy')
ims = ims[:,:38,:,:]
print(ims.shape)
for el in range(35641):
	ims[el,:,:,0] = gaussian_filter(ims[el,:,:,0],sigma=0.75)


ims = (ims-ims.min())/(ims.max()-ims.min())
mytr = ims[:33641,:,:,:]
myval = ims[33641:,:,:,:]
checkpointer = ModelCheckpoint(filepath='./leaky_edsr.h5',verbose=1,save_best_only=True)

history = mymodel.fit([data_tr],[mytr],batch_size=32,epochs=120,validation_data=([data_val],[myval]),callbacks=[checkpointer])
mymodel.save('leaky_edsr_final.h5')
my_loss = np.asarray(history.history['loss'])
my_val = np.asarray(history.history['val_loss'])
np.save('leaky_lr_echo_im_loss_no_code_hr.npy',my_loss)
np.save('leaky_lr_echo_im_val_no_code_hr.npy',my_val)




