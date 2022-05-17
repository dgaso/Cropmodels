from __future__ import print_function
import sys, os
from math import exp, log
from datetime import date

sys.path.append(r"C:\Users\gaso001\Thesis\campbell_diaz")
from pcse import signals
from pcse.base import SimulationObject, RatesTemplate, StatesTemplate, ParamTemplate, ParameterProvider
from pcse.traitlets import Float, Int, Instance, Enum,  Dict 
from pcse.decorators import prepare_rates, prepare_states
import pcse.exceptions as exc
from pcse.util import daylength, limit, AfgenTrait


class TemperatureReductionFactor(SimulationObject):

    alpha = Float()

    class Parameters(ParamTemplate):
        Topt = Float()
        Tmin = Float()
        Tmax = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        p = self.params
        self.alpha = log(2.)/(log((p.Tmax-p.Tmin)/(p.Topt-p.Tmin)))
        self._cache = {}

    def __call__(self, x, _cache={}):

        try:
            return _cache[x]
        except KeyError:
            pass

        p = self.params
        if x < p.Tmin or x > p.Tmax:
            v = _cache[x] = 0.
        else:
            p1 = (2*(x - p.Tmin)**self.alpha)
            p2 = ((p.Topt - p.Tmin)**self.alpha)
            p3 = ((x - p.Tmin)**(2*self.alpha))
            p4 = ((p.Topt - p.Tmin)**(2*self.alpha))
            v = _cache[x] = (p1 * p2 - p3)/p4

        return v


class PhotoperiodReductionFactor(SimulationObject):
    """Photoperiod reduction factor for soybean (short day)
       approach and parameters based on Setiyono et al. doi 10.1016/j.fcr.2006.07.011
       http://digitalcommons.unl.edu/agronomyfacpub/112
    """

    alpha = Float()
    m = Float(3)
    p0 = Float()

    class Parameters(ParamTemplate):
        MG = Float()  # Maturity Group rating
        Popt = Float()
        Pcrt = Float()

    def initialize(self, day, kiosk, parvalues):
        # First compute Popt and Pcrt based on maturity group rating
        # equation based on Setiyono et al. doi 10.1016/j.fcr.2006.07.011
        # http://digitalcommons.unl.edu/agronomyfacpub/112
        # if "Pcrt" in parvalues and "Popt" in parvalues:
            # # both optimal and critical daylength are there MG is not needed
            # # but dummy needs to be specified if missing (set to -99)
            # if "MG" not in parvalues:
                # parvalues["MG"] = -99
        # elif "MG" in parvalues:
            # MG = parvalues["MG"]
            # # derived Popt/Pcrt from MG based on setiyono equation
            # if "Popt" not in parvalues:
                # Popt = 12.759 - 0.388*MG - 0.058*MG**2
            # if "Pcrt" not in parvalues:
                # Pcrt = 27.275 - 0.493*MG - 0.066*MG**2
            # parvalues._cropdata.update({"Pcrt":Pcrt, "Popt":Popt})

        # derived Popt/Pcrt from MG based on setiyono equation
        #MG = parvalues["MG"]
        #Popt = 12.759 - 0.388*MG - 0.058*MG**2
        #Pcrt = 27.275 - 0.493*MG - 0.066*MG**2
        #parvalues.update({"Pcrt":Pcrt, "Popt":Popt})

        self.params = self.Parameters(parvalues)
        p = self.params
        self.alpha = log(2.)/log(((p.Pcrt - p.Popt)/self.m) + 1.)
        self.p0 = (p.Pcrt - p.Popt)/self.m

    def __call__(self, x):

        p = self.params
        if x < p.Popt:
            v =  1.
        elif x > p.Pcrt:
            v =  0.
        else:
            p1 = ((x - p.Popt)/self.m + 1.)
            p2 = ((p.Pcrt - x)/(p.Pcrt - p.Popt))
            v = (p1*(p2**self.p0))**self.alpha
        

        return v


class SoybeanPhenology(SimulationObject):
    """Implements the algorithms for phenologic development in WOFOST specifically for soybean.

    Phenologic development in WOFOST is expresses using a unitless scale which
    takes the values 0 at emergence, 1 at Anthesis (flowering) and 2 at
    maturity. This type of phenological development is mainly representative
    for cereal crops. All other crops that are simulated with WOFOST are
    forced into this scheme as well, although this may not be appropriate for
    all crops. For example, for potatoes development stage 1 represents the
    start of tuber formation rather than flowering.


    **Simulation parameters**

    =======  =============================================   =======  ============
     Name     Description                                     Type     Unit
    =======  =============================================   =======  ============
    TSUMEM   Temperature sum from sowing to emergence         SCr        |C| day
    TBASEM   Base temperature for emergence                   SCr        |C|
    TEFFMX   Maximum effective temperature for emergence      SCr        |C|
    DVRMAX1  Maximum development rate emergence to anthesis   SCr        |C| day
    DVRMAX2  Maximum develpment rate anthesis to maturity     SCr        |C| day
    DVSI     Initial development stage at emergence.          SCr        -
             Usually this is zero, but it can be higher
             for crops that are transplanted (e.g. paddy
             rice)
    DVSEND   Final development stage                          SCr        -
    MG       Maturity group rating for daylength sensivity    SCr        -
    Topt     Optimum temperature for phenological dev.        SCr        |C|
    Tmin     Temperature below which development is zero.     SCr        |C|
    Tmax     Temperature above which development is zero.     SCr        |C|
    =======  =============================================   =======  ============

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    DVS      Development stage                                  Y    -
    TSUM     Temperature sum                                    N    |C| day
    TSUME    Temperature sum for emergence                      N    |C| day
    DOS      Day of sowing                                      N    -
    DOE      Day of emergence                                   N    -
    DOR1     Day of R1 stage (beginning of flowering)           N    -
    DOR3     Day of R3 stage (pod development)                  N    -
    DOR5     Day of R5 stage (seed development)                 N    -
    DOR8     Day of R8 stage (fully ripe                        N    -
    STAGE    Current phenological stage, can take the           N    -
             folowing values:
             `emerging|vegetative|reproductive|mature`
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    DTSUME   Increase in temperature sum for emergence          N    |C|
    DTSUM    Increase in temperature sum for anthesis or        N    |C|
             maturity
    DVR      Development rate                                   Y    |day-1|
    =======  ================================================= ==== ============

    **External dependencies:**

    None

    **Signals sent or handled**

    `SoybeanPhenology` sends the `crop_finish` signal when maturity is
    reached and the `end_type` is 'maturity' or 'earliest'.

    """

    photoperiod_reduction_factor = Instance(PhotoperiodReductionFactor)
    temperature_reduction_factor = Instance(TemperatureReductionFactor)
    class Parameters(ParamTemplate):
        TSUMEM = Float(-99.)  # Temp. sum for emergence
        TBASEM = Float(-99.)  # Base temp. for emergence
        TEFFMX = Float(-99.)  # Max eff temperature for emergence
        DVRMAX1 = Float(-99.)  # Max development rate towards anthesis
        DVRMAX2 = Float(-99.)  # Max development rate towards maturity
        DVSI = Float(-99.)  # Initial development stage
        DVSEND = Float(-99.)  # Final development stage
        CROP_START_TYPE = Enum(["sowing", "emergence"])
        CROP_END_TYPE = Enum(["maturity", "harvest", "earliest"])

    #-------------------------------------------------------------------------------
    class RateVariables(RatesTemplate):
        DTSUME = Float(-99.)  # increase in temperature sum for emergence
        DVR = Float(-99.)  # development rate
        DAYL = Float(-99)
        PHOTORF = Float(-99)
        TEMPRF = Float(-99)

    #-------------------------------------------------------------------------------
    class StateVariables(StatesTemplate):
        DVS   = Float(-99.)  # Development stage
        TSUM  = Float(-99.)  # Temperature sum state
        TSUME = Float(-99.)  # Temperature sum for emergence state
        # States which register phenological events
        DOS = Instance(date)  # Day of sowing
        DOE = Instance(date)  # Day of emergence
        DOR1 = Instance(date)  # Day of start of flowering
        DOR3 = Instance(date)  # Day of pod development
        DOR5 = Instance(date)  # Day of seed filling
        DOR8 = Instance(date)  # Day of full ripeness
        DOH = Instance(date)  # Day of harvest
        STAGE = Enum([None, "emerging", "vegetative", "reproductive", "mature"])

    #---------------------------------------------------------------------------
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """
        self.photoperiod_reduction_factor = PhotoperiodReductionFactor(day, kiosk, parvalues)
        self.temperature_reduction_factor = TemperatureReductionFactor(day, kiosk, parvalues)
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk

        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)

        # Define initial states
        DVS = self.params.DVSI
        DOS, DOE, STAGE = self._get_initial_stage(day)
        self.states = self.StateVariables(kiosk, publish="DVS",
                                          TSUM=0., TSUME=0., DVS=DVS,
                                          DOS=DOS, DOE=DOE, DOR1=None,
                                          DOR3=None, DOR5=None, DOR8=None,
                                          DOH=None, STAGE=STAGE)

    def _get_initial_stage(self, day):
        """"""
        p = self.params

        # Define initial stage type (emergence/sowing) and fill the
        # respective day of sowing/emergence (DOS/DOE)
        if p.CROP_START_TYPE == "emergence":
            STAGE = "vegetative"
            DOE = day
            DOS = None

            # send signal to indicate crop emergence
            self._send_signal(signals.crop_emerged)

        elif p.CROP_START_TYPE == "sowing":
            STAGE = "emerging"
            DOS = day
            DOE = None

        else:
            msg = "Unknown start type: %s" % p.CROP_START_TYPE
            raise exc.PCSEError(msg)

        return DOS, DOE, STAGE

    @prepare_rates
    def calc_rates(self, day, drv):
        """Calculates the rates for phenological development
        """
        p = self.params
        r = self.rates
        s = self.states

        # Day length
        r.DAYL = daylength(day, drv.LAT, angle=-0.83)

        # temperature reduction factor
        r.TEMPRF = self.temperature_reduction_factor(drv.TEMP)

        # photoperiod
        r.PHOTORF = 1.0

        if s.STAGE == "emerging":
            r.DTSUME = limit(0., (p.TEFFMX - p.TBASEM), (drv.TEMP - p.TBASEM))
        elif s.STAGE == 'vegetative':
            r.DVR = p.DVRMAX1 * r.TEMPRF
        elif s.STAGE in ['reproductive', 'mature']:
            # photoperiod reduction factor
            r.PHOTORF = self.photoperiod_reduction_factor(r.DAYL)
            r.DVR = p.DVRMAX2 * r.PHOTORF * r.TEMPRF
        else: # Problem: no stage defined
            msg = "No STAGE defined in phenology submodule"
            raise exc.PCSEError(msg)

        msg = "Finished rate calculation for %s"
        self.logger.debug(msg % day)

    #---------------------------------------------------------------------------
    @prepare_states
    def integrate(self, day, delt):
        """Updates the state variable and checks for phenologic stages
        """

        p = self.params
        r = self.rates
        s = self.states

        if s.STAGE == "emerging":
            s.TSUME += r.DTSUME
            if s.TSUME >= p.TSUMEM:
                self._next_stage(day)

        elif s.STAGE == 'vegetative':
            s.DVS += r.DVR
            if s.DVS >= 1.0:
                self._next_stage(day)
                s.DVS = 1.0

        elif s.STAGE == 'reproductive':
            s.DVS += r.DVR
            # Check of R5 stage is reached at DVS=1.15.
            if s.DVS >= 1.15 and s.DOR5 is None:
                s.DOR5 = day

            # Check if maturity is reached
            if s.DVS >= p.DVSEND:
                self._next_stage(day)

        elif s.STAGE == 'mature':
            s.DVS += r.DVR


        else: # Problem no stage defined
            msg = "No STAGE defined in phenology submodule"
            raise exc.PCSEError(msg)

        msg = "Finished state integration for %s"
        self.logger.debug(msg % day)

    #---------------------------------------------------------------------------
    def _next_stage(self, day):
        """Moves states.STAGE to the next phenological stage"""
        s = self.states
        p = self.params

        current_STAGE = s.STAGE
        if s.STAGE == "emerging":
            s.STAGE = "vegetative"
            s.DOE = day
            # send signal to indicate crop emergence
            self._send_signal(signals.crop_emerged)

        elif s.STAGE == "vegetative":
            s.STAGE = "reproductive"
            s.DOR1 = day

        elif s.STAGE == "reproductive":
            s.STAGE = "mature"
            s.DOR8 = day
            if p.CROP_END_TYPE in ["maturity","earliest"]:
                self._send_signal(signal=signals.crop_finish,
                                  day=day, finish_type="maturity")
        elif s.STAGE == "mature":
            msg = "Cannot move to next phenology stage: maturity already reached!"
            raise exc.PCSEError(msg)

        else: # Problem no stage defined
            msg = "No STAGE defined in phenology submodule."
            raise exc.PCSEError(msg)

        msg = "Changed phenological stage '%s' to '%s' on %s"
        self.logger.info(msg % (current_STAGE, s.STAGE, day))

    #---------------------------------------------------------------------------
    def _on_CROP_FINISH(self, day, finish_type=None):
        """Handler for setting day of harvest (DOH). Although DOH is not
        strictly related to phenology (but to management) this is the most
        logical place to put it.
        """
        if finish_type == 'harvest':
            self._for_finalize["DOH"] = day
