from campbell_diaz.model import CampbellDiazModel
import numpy as np
import pandas as pd

class ModelRerunner(object):
    """Reruns a given model with different values of parameters rdmax, fc and pwp.

    Returns a pandas DataFrame with simulation results of the model with given
    parameter values.
    """
    parameters = ["WUE", "RDMAX","FNTR","initLAI"]# ,

    def __init__(self, params, wdp, agro):
        self.params = params
        self.wdp = wdp
        self.agro = agro

    def __call__(self, par_values):
        # Check if correct number of parameter values were provided
        if len(par_values) != len(self.parameters):
            msg = "Optimizing %i parameters, but only % values were provided!" % (len(self.parameters, len(par_values)))
            raise RuntimeError(msg)
        # Clear any existing overrides
        self.params.clear_override()
        # Set overrides for the new parameter values
        for parname, value in zip(self.parameters, par_values):
            self.params.set_override(parname, value)
        #self.params.set_override(varname='FCP', value=0.35,check=True)
        #self.params.set_override(varname='PWPP', value=0.33,check=True)  

        # Run the model with given parameter values
        engine = CampbellDiazModel(self.params, self.wdp, self.agro)
        engine.run_till_terminate()
        df = pd.DataFrame(engine.get_output())
        df.index = pd.to_datetime(df.day)
        return df


##Define an object function calculator
class ObjectiveFunctionCalculator(object):
    """Computes the objective function.

    This class runs the simulation model with given parameter values and returns the objective
    function as the sum of squared difference between observed and simulated LAI.
.   """

    def __init__(self, params, wdp, agro, observations):
        self.modelrerunner = ModelRerunner(params, wdp, agro)
        self.df_observations = observations
        self.n_calls = 0

    def __call__(self, par_values, grad=None):
        """Runs the model and computes the objective function for given par_values.

        The input parameter 'grad' must be defined in the function call, but is only
        required for optimization methods where analytical gradients can be computed.
        """
        self.n_calls += 1
        print(".", end="")
        # Run the model and collect output
        self.df_simulations = self.modelrerunner(par_values)
        # compute the differences by subtracting the DataFrames
        # Note that the dataframes automatically join on the index (dates) and column names
        combined = self.df_observations.join(self.df_simulations, how="left", rsuffix="_sim")
        df_differences = combined.LAI - combined.LAI_sim
        # Compute the RMSE on the LAI column
        obj_func = np.sqrt(np.mean(df_differences ** 2))
        self.err=obj_func
        rel_error = (self.err/np.mean(combined.LAI))*100
        self.rrmse = rel_error
        LAImax = np.max(combined.LAI)
        self.LAImx = LAImax

        return obj_func


