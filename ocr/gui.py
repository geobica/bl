import seabreeze
from matplotlib import pyplot as plt
import numpy as np
import time
import pickle
from scipy.constants import h, k
from scipy.constants import c as c_light
from seabreeze.spectrometers import Spectrometer
import os

import multiprocessing as mp
import sys
import numpy as np
import random

from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QGridLayout,
    QFrame,
    QHBoxLayout,
    QVBoxLayout,
    QCheckBox,
    QLabel,
    QLineEdit,
)
from PyQt5.QtCore import QTimer
from PyQt5.QtGui import QColor

import pyqtgraph as pg


THROUGH_CIEXYZ, WORKING_MODEL = 0, 1

# the method to be used for reconstructing RGB colors from a spectrum
conversion_method_used = THROUGH_CIEXYZ


# this was recorded with the cap on, it should be all 0, but it isn't, so this can be used to calibrate.
CALIB_with_cap_on_S = np.load("CALIB_with_cap_on_S.npy")
recorded_calibration_baseline = np.mean(CALIB_with_cap_on_S,axis=0)

# this was a recording of a white sheet of paper illuminated outside under cloud cover
# to do a real calibration I would need to get a recording of a true blackbody, or at least a recording of an more even reflector under the sun
# this works for my basic purposes for now, but it is not a scientifically rigorous calibration
CALIB_white_sheet_S = np.load("CALIB_white_sheet_S.npy")
recorded_calibration_baseline = np.mean(CALIB_with_cap_on_S,axis=0)
white_sheet_baseline = np.min(np.mean(CALIB_white_sheet_S[:,:],axis=0)-recorded_calibration_baseline)
recorded_calibration_whites = np.mean(CALIB_white_sheet_S[:,:],axis=0)-recorded_calibration_baseline-white_sheet_baseline
# to make sure this can be used for division without values being abnormally high
recorded_calibration_whites = np.maximum(100,recorded_calibration_whites)


def planck(lam, T):
    lam_m = lam / 1.0e9
    fac = h * c_light / lam_m / k / T
    B = 2 * h * c_light**2 / lam_m**5 / (np.exp(fac) - 1)
    return B


c = time.time()
seabreeze.use("pyseabreeze")


spec = Spectrometer.from_first_available()
# the first ten are skipped because they are unusually low
wavelengths = spec.wavelengths()[10:]

lam = spec.wavelengths()[10:]

blackbody = planck(lam, 6500)
blackbody /= np.max(blackbody)
# plt.plot(lam,blackbody,label="6500K")
# plt.plot(lam,recorded_calibration_whites/np.max(recorded_calibration_whites),label="white sheet")
# plt.legend()
# plt.show()


def probe():
    integration_time = 20000 # default
    try:
        with open("integration_time.txt","r") as r:
            # outside of these bounds seabreeze breaks
            integration_time = max(3000,min(655350000,int(r.read())))
    except:
        pass
    spec.integration_time_micros(integration_time)
    intensity_total = spec.intensities().astype(float)
    return spec.wavelengths()[10:], intensity_total[10:]

# for THROUGH_CIEXYZ
# this file contains CIE 1931 XYZ color matching functions
# https://www.site.uottawa.ca/~edubois/mdsp/data/ciexyz31_1.txt
import pandas as pd
ciexyz31_1 = pd.read_csv("ciexyz31_1.txt",header=None,names=["λ","X","Y","Z"])
ciexyz_data = np.array(ciexyz31_1.values)

ciexyz_data_lam = np.zeros((lam.shape[0], 4))
for i in range(lam.shape[0]):
    try:
        ciexyz_data_lam[i] = ciexyz_data[np.where(ciexyz_data[:, 0] < lam[i])[0][-1]]
    except:
        pass
        
# for WORKING_MODEL
working_model = None
working_model_E = None
if os.path.exists("models/working_model_1767374906404002560.sav"):
    working_model = pickle.load(open("models/working_model_1767374906404002560.sav", "rb"))
    working_model_E = np.load("models/working_model_1767374906404002560_E.npy")


def spectrum_to_rgb(data):
    if conversion_method_used==THROUGH_CIEXYZ:
        rad_block = np.maximum(0, np.array(data))
        xyz_sums = np.sum(rad_block[:, None] * ciexyz_data_lam[:, 1:], axis=0)

        xyz_srgb = np.array(
            [
                [3.2404542, -1.5371385, -0.4985314],
                [-0.9692660, 1.8760108, 0.041556],
                [0.0556434, -0.2040259, 1.0572252],
            ]
        )
        srgb = np.matmul(xyz_srgb, xyz_sums[:, None])[
            :, 0
        ]
        multipliers = np.exp(np.linspace(-2, 2, num=9, endpoint=True))
        srgb = srgb[:, None] * multipliers[None, :]
        below_threshold = srgb < 0.0031308
        srgb[below_threshold] = 12.92 * srgb[below_threshold]
        srgb[np.logical_not(below_threshold)] = (
            1.055 * srgb[np.logical_not(below_threshold)] ** (1 / 2.4) - 0.055
        )
        new_three_channel = np.maximum(
            0, np.minimum(255, srgb / np.max(srgb) * 255)
        ).astype(int)
        return new_three_channel.T
    if conversion_method_used==WORKING_MODEL:
        multipliers = np.exp(np.linspace(-2, 2, num=9, endpoint=True))
        v = np.matmul(working_model_E, data[:, None] * multipliers[None, :])[:3].T
        return np.maximum(0, np.minimum(255, working_model.predict(v))).astype("uint8")




# Background probe process
def probe_worker(queue: mp.Queue):
    """
    Continuously run probe() and push results to the queue.
    """
    while True:
        wav, data = probe()
        if data is not None:
            queue.put(data)




class ColorWidget(QFrame):
    """
    Widget to display a randomly generated color for calibration
    """
    def __init__(self):
        super().__init__()
        self.setFrameStyle(QFrame.Panel | QFrame.Raised)
        self.setAutoFillBackground(True)

    def set_random_color(self):
        self.this_parent.color_index = (self.this_parent.color_index + 619) % 1000

        x = []
        y = []
        failed = False
        try:
            color_index = self.this_parent.history[-1, 0]

            where_this_color = np.where(self.this_parent.history[:, 0] == color_index)[
                0
            ]
            where_to_use = where_this_color[int(len(where_this_color) / 2) :]
            if len(where_to_use) > 3:
                if (
                    int(color_index) in self.this_parent.color_record
                    and int(color_index) > 0
                ):
                    x.append(np.array(self.this_parent.color_record[color_index]))
                    y.append(
                        np.mean(self.this_parent.history[where_to_use], axis=0)[10:]
                    )
            if not len(x):
                failed = True
        except:
            failed = True
        if not failed:
            # store the spectrum mean
            self.this_parent.S.append(y[-1])

            # store the RGB of the color
            self.this_parent.T.append(x[-1])

        rgb = [random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)]
        color = QColor(rgb[0], rgb[1], rgb[2])
        self.this_parent.color_record[int(self.this_parent.color_index)] = rgb
        palette = self.palette()
        palette.setColor(self.backgroundRole(), color)
        self.setPalette(palette)


class OtherColorWidget(QFrame):
    """
    Widget to display a reconstructed color from what the spectrometer is seeing
    """
    def __init__(self):
        super().__init__()
        self.setFrameStyle(QFrame.Panel | QFrame.Raised)
        self.setAutoFillBackground(True)

        self.image_plot = pg.PlotWidget(title="Reconstructed Colors")
        self.image_item = pg.ImageItem()
        self.image_plot.addItem(self.image_item)
        self.image_plot.setMouseEnabled(x=False, y=False)

        layout = QGridLayout(self)
        layout.setSpacing(5)
        layout.addWidget(self.image_plot, 0, 0)

    def set_random_color(self):
        # spectrum_to_rgb returns the color at 9 different multiplier values to simulate different brightness levels
        usable_history = self.this_parent.history[
            np.mean(self.this_parent.history, axis=1) > 200
        ]
        if len(usable_history) > 10:
            baseline_removed = np.mean(
                self.this_parent.history[-5:, 10:].astype(float), axis=0
            )-recorded_calibration_baseline - np.percentile(
                usable_history[: max(1, len(usable_history) - 20), 10:].astype(float)-recorded_calibration_baseline,
                10,
            )
            the_colors = spectrum_to_rgb(blackbody * baseline_removed/recorded_calibration_whites)

            self.image_item.setImage(the_colors.reshape((3, 3, 3)))

            the_color = the_colors[4]
            if the_color[0] >= 0 and the_color[0] <= 255:
                color = QColor(
                    the_color[0],
                    the_color[1],
                    the_color[2],
                )
                palette = self.palette()
                palette.setColor(self.backgroundRole(), color)
                self.setPalette(palette)


class MainWindow(QWidget):
    def __init__(self, data_queue: mp.Queue):
        super().__init__()

        self.queue = data_queue
        self.data = np.zeros((2038,))

        self.T = []  # RGB colors
        self.S = []  # spectrums seen when looking at those colors displayed

        self.color_index = 0
        self.color_record = {}

        self.setWindowTitle("The Monitor")
        self.resize(1200, 700)

        # Main vertical layout
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(8)

        # Layout of bar at top
        top_layout = QHBoxLayout(self)
        top_layout.setSpacing(8)

        # Title
        title_label = QLabel("OceanOptics Spectrometer Monitor")
        title_label.setStyleSheet("font-weight: bold; font-size: 14px;")
        top_layout.addWidget(title_label)

        input_label = QLabel("Integration Time (μs):")
        top_layout.addWidget(input_label)

        # Text input above the grid
        self.integration_time = 20000 # default
        self.collect_rate_input = QLineEdit()
        self.collect_rate_input.setPlaceholderText("20000")
        self.collect_rate_input.textChanged.connect(self.change_integration_time)
        top_layout.addWidget(self.collect_rate_input)

        # original grid layout
        grid_layout = QGridLayout()
        grid_layout.setSpacing(5)

        main_layout.addLayout(top_layout)
        main_layout.addLayout(grid_layout)

        # Live line plot (top-left)
        self.line_plot = pg.PlotWidget(title="Live Spectrum")
        self.line_curve = self.line_plot.plot(wavelengths, wavelengths, pen="y")
        self.line_curve_mean = self.line_plot.plot(
            wavelengths, np.ones_like(wavelengths) * 0, pen="w"
        )
        # self.line_curve_R = self.line_plot.plot(
        #     wavelengths, np.ones_like(wavelengths) * 400, pen="r"
        # )
        # self.line_curve_G = self.line_plot.plot(
        #     wavelengths, np.ones_like(wavelengths) * 400, pen="g"
        # )
        # self.line_curve_B = self.line_plot.plot(
        #     wavelengths, np.ones_like(wavelengths) * 400, pen="b"
        # )
        # self.line_curve_k = self.line_plot.plot(
        #     wavelengths, np.ones_like(wavelengths) * 400, pen="k"
        # )
        # self.line_curve_w = self.line_plot.plot(
        #     wavelengths, np.ones_like(wavelengths) * 400, pen="w"
        # )

        # Scrolling image on right to show recorded data
        self.image_plot = pg.PlotWidget(title="Recording")
        self.image_item = pg.ImageItem()
        self.image_plot.addItem(self.image_item)
        self.image_plot.setMouseEnabled(x=False, y=False)

        self.history_length = 800
        self.data_length = 2048
        self.history = np.zeros((self.history_length, self.data_length))
        self.history_indices = np.zeros((self.history_length,))
        self.history_timestamps = np.zeros((self.history_length,))

        # The random and reconstructed colors
        self.color_widget = ColorWidget()
        self.color_widget.this_parent = self
        self.other_color_widget = OtherColorWidget()
        self.other_color_widget.this_parent = self

        # Grid widgets
        grid_layout.addWidget(self.line_plot, 0, 0, 1, 2)
        grid_layout.addWidget(self.image_plot, 0, 2, 2, 2)
        grid_layout.addWidget(self.color_widget, 1, 0)
        grid_layout.addWidget(self.other_color_widget, 1, 1)

        grid_layout.setColumnStretch(0, 1)
        grid_layout.setColumnStretch(1, 2)
        grid_layout.setRowStretch(0, 1)
        grid_layout.setRowStretch(1, 1)

        # Timers
        self.gui_timer = QTimer()
        self.gui_timer.timeout.connect(self.consume_queue)
        self.gui_timer.start(30) # fast GUI update loop

        self.color_timer = QTimer()
        self.color_timer.timeout.connect(self.color_widget.set_random_color)
        self.color_timer.start(3000) # get a new random color every 3 seconds

        self.other_color_timer = QTimer()
        self.other_color_timer.timeout.connect(self.other_color_widget.set_random_color)
        self.other_color_timer.start(30)

    def change_integration_time(self):
        try:
            self.integration_time = int(self.collect_rate_input.text());
        except Exception as e:
            print(e)
            self.integration_time = 20000
        with open("integration_time.txt","w") as f:
            f.write(f"{self.integration_time}")
        print("self.integration_time",self.integration_time)

    def consume_queue(self):
        """
        Drain the queue and display the most recent data.
        """
        latest = None
        while not self.queue.empty():
            latest = self.queue.get_nowait()

        if latest is None:
            return

        data = np.asarray(latest, dtype=float)
        self.data = np.array(data)

        # Update line plot
        self.line_curve.setData(wavelengths, data - np.percentile(
                data.astype(float),
                10,
            ))
        usable_history = self.history[
            np.mean(self.history, axis=1) > 200
        ]
        if len(usable_history) > 10:
            baseline_removed = np.mean(
                self.history[-5:, 10:].astype(float), axis=0
            )-recorded_calibration_baseline - np.percentile(
                usable_history[: max(1, len(usable_history) - 20), 10:].astype(float)-recorded_calibration_baseline,
                10,
            )
            self.line_curve_mean.setData(wavelengths, 1000*blackbody/recorded_calibration_whites * baseline_removed)

        # Update scrolling image
        self.history = np.roll(self.history, -1, axis=0)
        self.history_indices = np.roll(self.history_indices, -1, axis=0)
        self.history_timestamps = np.roll(self.history_timestamps, -1, axis=0)
        self.history[-1, 10:] = data
        self.history_indices[-1] = self.history_indices[-2] + 1
        self.history_timestamps[-1] = time.time_ns()
        if self.history_indices[-1] % self.history_length == self.history_length - 1:
            all_history = np.array(self.history).astype("int64")
            all_history[:, 0] = self.history_indices
            all_history[:, 1] = self.history_timestamps
            np.save(f"recordings/{int(self.history_timestamps[-1])}.npy", all_history)

        # The bottom 10 are used to show the randomly generated colors
        self.history[-1, :10] = self.color_index
        rgb_history = np.zeros((self.history.shape[0], self.history.shape[1], 3))
        for i in range(3):
            rgb_history[:, :, i] = np.maximum(
                0, np.minimum(255, (self.history / 1000 * 255))
            )
        for color_key in np.unique(self.history[:, 0]):
            where_this_color = np.where(self.history[:, 0] == color_key)[0]
            if int(color_key) in self.color_record and int(color_key) > 0:
                the_color = self.color_record[int(color_key)]
                for i in range(3):
                    history_where = rgb_history[where_this_color]
                    history_where[:, :10, i] = the_color[i]
                    rgb_history[where_this_color] = history_where

        self.image_item.setImage(rgb_history.astype("uint8"))

    def closeEvent(self, *args, **kwargs):
        print("Closing!!")
        all_history = np.array(self.history).astype("int64")
        all_history[:, 0] = self.history_indices
        all_history[:, 1] = self.history_timestamps
        # record important data
        np.save(
            f"recordings/{int(self.history_timestamps[-1])}_T.npy", np.array(self.T)
        )
        np.save(
            f"recordings/{int(self.history_timestamps[-1])}_S.npy", np.array(self.S)
        )
        np.save(f"recordings/{int(self.history_timestamps[-1])}_c.npy", all_history)


if __name__ == "__main__":
    mp.set_start_method("spawn")

    data_queue = mp.Queue(maxsize=5)

    probe_process = mp.Process(target=probe_worker, args=(data_queue,), daemon=True)
    probe_process.start()

    app = QApplication(sys.argv)
    window = MainWindow(data_queue)
    window.show()

    exit_code = app.exec_()
    probe_process.terminate()
    sys.exit(exit_code)
