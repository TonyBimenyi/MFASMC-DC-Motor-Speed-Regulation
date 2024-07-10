import numpy as np

from Controllers.MFASM_controller import MFASM_controller
from Reference.trajectory import generate_trajectory
from Output.output import Output
from Estimator.estimator import initialize_phi
from Plot.plotter import plot_results