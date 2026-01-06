import os

import numpy as np
from matplotlib import pyplot as plt
import scipy
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn import model_selection
import pickle

most_recent = sorted(os.listdir("recordings"))[-3:]

S = np.load(f"recordings/{most_recent[0]}")
T = np.load(f"recordings/{most_recent[1]}")

if False:
    # for some reason in an earlier run there was a trend where it started with a higher baseline and drifted lower by about 11 DN over the hour and a half it was running
    # this code is to show that drift
    plt.title("Baseline Drift")
    plt.plot(np.mean(S, axis=1), label="mean of each signal")
    plt.plot(scipy.signal.medfilt(np.mean(S, axis=1), 9), label="medfilt of means")
    plt.legend()
    plt.show()

baselines = []
for i in range(len(S)):
    baselines.append(
        np.percentile(np.mean(S[max(i - 50, 0) : min(i + 50, len(S) - 1)], axis=0), 10)
    )

# remove the drifting baseline
S = S - np.array(baselines)[:, None]

if False:
    # plot to show that both the highs and lows are essentially random after removing the drifting baseline
    lows = []
    highs = []
    for i in range(int(len(S) / 100)):
        lows.append(np.percentile(np.mean(S[100 * i : 100 * i + 100], axis=0), 10))
        highs.append(np.percentile(np.mean(S[100 * i : 100 * i + 100], axis=0), 90))

    plt.plot(lows, highs)
    plt.scatter(lows, highs)
    plt.show()

x = np.hstack((np.array(T), np.ones((np.array(T).shape[0]))[:, None]))
y = np.array(S)
# x is nx4
# y is nx2038
# where n is the number of colors recorded
# E is 4x2038
# y = x E

# calculate the best 1000x4 array to use to reconstruct
lstsq_res = np.linalg.lstsq(x, y, rcond=None)
print(lstsq_res)
E = lstsq_res[0]
plt.plot(E[0], label="Red")
plt.plot(E[1], label="Green")
plt.plot(E[2], label="Blue")
plt.plot(E[3], label="Baseline")
plt.legend()
plt.show()


# but the real interesting part is getting y from x
# ie getting the RGB values for a spectrum
# convolution is v = E S
print(E.shape, S.shape)
v = np.matmul(E, S.T)[:3].T
print("v", v.shape)
print("x", x.shape)
lstsq_res_M = np.linalg.lstsq(v, T, rcond=None)
print(lstsq_res_M)
M = lstsq_res_M[0]

C = np.matmul(v, M)
rgb_im = np.dstack((C, T)).transpose((0, 2, 1))
rgb_im = np.repeat(rgb_im, 100, axis=1)
plt.imshow(np.maximum(0, np.minimum(255, rgb_im)).astype("uint8"))
plt.show()

# Build pipeline: scaling + nonlinear model
model = Pipeline(
    [
        ("scaler", StandardScaler()),
        (
            "mlp",
            MLPRegressor(
                hidden_layer_sizes=(256, 16),
                activation="relu",
                solver="adam",
                max_iter=500,
                random_state=0,
            ),
        ),
    ]
)
seed = 7
X_train, X_test, Y_train, Y_test = model_selection.train_test_split(
    v, T, test_size=0.2, random_state=seed
)
# Fit nonlinear inverse
model.fit(X_train, Y_train)
# save the model to disk
filename = f'models/working_model_{most_recent[0].split("_")[0]}.sav'
pickle.dump(model, open(filename, "wb"))
np.save(f'models/working_model_{most_recent[0].split("_")[0]}_E.npy', E)

C = model.predict(v)
rgb_im = np.dstack((C, T)).transpose((0, 2, 1))
rgb_im = np.repeat(rgb_im, 100, axis=1)
plt.title("distance from truth to predicted")
plt.plot(np.linalg.norm(C - T, axis=1), label="")
plt.show()
plt.imshow(np.maximum(0, np.minimum(255, rgb_im)).astype("uint8"))
plt.show()


print(M, M.shape)

# plt.plot
for recording_file in os.listdir("recordings"):
    recording_arr = np.load(f"recordings/{recording_file}")
    plt.plot(recording_arr[:, 0])
    plt.show()
