
#!/usr/bin/env python

import numpy as np
from tensorflow import keras
from tensorflow.keras import layers
import pandas as pd

"""
## Prepare the data
"""

# Model / data parameters
num_classes = 2#10
#input_shape = (28,28,1)
input_shape = (574, 499, 1)

# the data, split between train and test sets
#(x_train, y_train), (x_test, y_test) = keras.datasets.mnist.load_data()
#print("x_train shape:", x_train.shape)
#print("y_train shape:", y_train.shape)
#print("x_test shape:", x_test.shape)
#print("y_test shape:", y_test.shape)

# read four XLS files: SamplesTr, SelectionTr, SamplesTe, SelectionTe
x_train = pd.read_csv('~/genomedk/DELFI2/Workspaces/matov/MNIST/samplesTr.csv') # 69+64 x 574 x 499
y_train = pd.read_csv('~/genomedk/DELFI2/Workspaces/matov/MNIST/selectionTr.csv') # vector 69 1s and 64 0s
x_test = pd.read_csv('~/genomedk/DELFI2/Workspaces/matov/MNIST/samplesTe.csv') # 10+10 x 574 x 499
y_test = pd.read_csv('~/genomedk/DELFI2/Workspaces/matov/MNIST/selectionTe.csv') # vector 10 1s and 10 0s
#x_train = pd.read_csv('~/matovanalysis/DELFI_analysis/python/samplesTr.csv') # 69+64 x 574 x 499
#y_train = pd.read_csv('~/matovanalysis/DELFI_analysis/python/selectionTr.csv') # vector 69 1s and 64 0s
#x_test = pd.read_csv('~/matovanalysis/DELFI_analysis/python/samplesTe.csv') # 10+10 x 574 x 499
#y_test = pd.read_csv('~/matovanalysis/DELFI_analysis/python/selectionTe.csv') # vector 10 1s and 10 0s
print("x_train shape:", x_train.shape)
print("y_train shape:", y_train.shape)
print("x_test shape:", x_test.shape)
print("y_test shape:", y_test.shape)
# Scale images to the [0, 1] range
#x_train = x_train.astype("float32") / 255
#x_test = x_test.astype("float32") / 255
# Make sure images have shape (28, 28, 1)
#x_train = np.expand_dims(x_train, -1)
#x_test = np.expand_dims(x_test, -1)
#print("x_train shape:", x_train.shape)
#print(x_train.shape[0], "train samples")
#print(x_test.shape[0], "test samples")


# convert class vectors to binary class matrices
#y_train = keras.utils.to_categorical(y_train, num_classes)
#y_test = keras.utils.to_categorical(y_test, num_classes)

"""
## Build the model
"""

model = keras.Sequential(
    [
        keras.Input(shape=input_shape),
        layers.Conv2D(32, kernel_size=(3, 3), activation="relu"),
        layers.MaxPooling2D(pool_size=(2, 2)),
        layers.Conv2D(64, kernel_size=(3, 3), activation="relu"),
        layers.MaxPooling2D(pool_size=(2, 2)),
        layers.Flatten(),
        layers.Dropout(0.5),
        layers.Dense(num_classes, activation="softmax"),
    ]
)

model.summary()

"""
## Train the model
"""

batch_size = 128
epochs = 15

model.compile(loss="categorical_crossentropy", optimizer="adam", metrics=["accuracy"])

model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs, validation_split=0.1)

"""
## Evaluate the trained model
"""

score = model.evaluate(x_test, y_test, verbose=0)
print("Test loss:", score[0])
print("Test accuracy:", score[1])
