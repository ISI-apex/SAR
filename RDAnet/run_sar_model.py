import tensorflow as tf
import numpy as np
import sys

mymodel = tf.keras.models.load_model(sys.argv[1]+'.h5')

mymodel.summary()


ims = np.zeros((35641,256,256,2),dtype=np.float32)
ims[:33641,:,:,:] = np.load('/home/echo_hr_108_train.npy')
ims[33641:,:,:,:] = np.load('/home/echo_hr_108_val.npy')
data_mean = np.load('/home/echo_hr_108_mean.npy')

ims = (ims-data_mean)

dout = mymodel.predict(ims)

np.save(sys.argv[1]+'_results.npy',dout)
