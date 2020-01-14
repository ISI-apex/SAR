import time
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import GaussianDropout, Dense, Permute, Dropout, Activation, Reshape, Flatten, GaussianNoise
from tensorflow.keras.layers import Input, SpatialDropout2D, Conv2D, Conv2DTranspose, BatchNormalization, AveragePooling2D, MaxPooling2D
from tensorflow.keras.optimizers import SGD, Adam, Adadelta
from tensorflow.keras import regularizers
from tensorflow.keras import backend as K
from tensorflow.keras.models import Model
import h5py
import numpy as np
from copy import deepcopy

c0 = Input(shape=(None,None,1)) 
c1 = Conv2D(filters=32, kernel_size=(3,3), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv1')(c0)
c2 = Conv2D(filters=32, kernel_size=(3,3), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv2')(c1)
c3 = MaxPooling2D(pool_size=(2,2),padding='same',data_format='channels_last')(c2)

c4 = Conv2D(filters=64, kernel_size=(3,3), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv3')(c3)
c5 = Conv2D(filters=64, kernel_size=(3,3), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv4')(c4)
c6 = MaxPooling2D(pool_size=(2,2),padding='same',data_format='channels_last')(c5)

c7 = Conv2D(filters=128, kernel_size=(3,3), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv5')(c6)
c8 = Conv2D(filters=128, kernel_size=(3,3), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv6')(c7)
c9 = MaxPooling2D(pool_size=(2,2),data_format='channels_last',padding='same')(c8)

c10 = Conv2D(filters=256, kernel_size=(3,3), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv7')(c9)
c11 = Conv2D(filters=256, kernel_size=(3,3), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv8')(c10)
c12 = MaxPooling2D(pool_size=(2,2),data_format='channels_last',padding='same')(c11)

c13 = Conv2D(filters=512, kernel_size=(4,4), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation='relu',activity_regularizer=None,kernel_initializer='glorot_normal',name='conv9')(c12)

#c14 = GaussianDropout(0.5)(c13)
c14 = SpatialDropout2D(0.5)(c13)
c15 = Conv2D(filters=11, kernel_size=(1,1), strides=(1,1),use_bias=True,  data_format="channels_last",padding="same",activation=None,activity_regularizer=None,kernel_initializer='glorot_normal',name='conv10')(c14)

c16 = Conv2DTranspose(filters=11, kernel_size=(32,32), strides=(16,16), data_format="channels_last",padding="same",activation=None,kernel_initializer='glorot_normal',name='conv11' )(c15)

c17 = Reshape((-1,11))(c16)
c18 = Activation('softmax')(c17)

rms = Adam(lr=1e-5)
mymodel = Model(inputs=c0,outputs=c18)
mymodel.compile(loss='sparse_categorical_crossentropy', optimizer=rms,metrics=['sparse_categorical_accuracy'])


print('ready to goo!!!')
mymodel.summary()
mymodel.load_weights('sar_atr.h5')




myf = np.load('my_192x192f.npy')

yout = mymodel.predict(myf,batch_size=1)
yout2 = yout[0,:,:]
yout3 = np.argmax(yout2,axis=1)
yout4 = yout3.reshape((192,192))

np.save('atr_im.npy',yout4)








