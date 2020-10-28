import ctypes
import numpy as np
import matplotlib.pyplot as plt

class RPASim():
    '''docstring for RPASim.'''

    def __init__(self, path):
        # Load wrapper library
        self.rpa = ctypes.CDLL(f'{path}/libwrapper.dll')

        # # Declare used functions
        # self.rpa.initialize.argtypes = [ctypes.c_bool]
        self.rpa.initializeWithPath.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_bool]
        #
        self.rpa.configFile.restype = ctypes.c_void_p
        self.rpa.configFileSaveAs.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        self.rpa.configFileCombustionChamberConditionsSetPressure.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFileNozzleFlowOptionsSetCalculateNozzleFlow.argtypes = [ctypes.c_void_p, ctypes.c_bool]
        self.rpa.configFileNozzleFlowOptionsSetNozzleInletConditions.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFileNozzleFlowOptionsSetNozzleExitConditions.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFilePropellantAddOxidizer.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFilePropellantAddFuel.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFilePropellantSetRatio.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFileEngineSizeSetThrust.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFileEngineSizeSetAmbientPressure.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFileEngineSizeSetMdot.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFileEngineSizeSetThroatDiameter.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.configFileEngineSizeSetChamberLength.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p, ctypes.c_int]
        self.rpa.configFileEngineSizeSetContractionAngle.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p]

        self.rpa.performanceCreate.argtypes = [ctypes.c_void_p, ctypes.c_bool, ctypes.c_bool]
        self.rpa.performanceCreate.restype = ctypes.c_void_p
        self.rpa.performanceDelete.argtypes = [ctypes.c_void_p]
        self.rpa.performanceSolve.argtypes = [ctypes.c_void_p, ctypes.c_bool]
        self.rpa.performanceGetOF.argtypes = [ctypes.c_void_p]
        self.rpa.performanceGetOF.restype = ctypes.c_double
        self.rpa.performanceGetDeliveredIsp.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_char_p, ctypes.c_double]
        self.rpa.performanceGetDeliveredIsp.restype = ctypes.c_double
        self.rpa.performanceGetDeliveredIspH.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_char_p, ctypes.c_double]
        self.rpa.performanceGetDeliveredIspH.restype = ctypes.c_double
        self.rpa.performanceGetIdealIsp.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.performanceGetIdealIsp.restype = ctypes.c_double
        self.rpa.performanceGetIdealIspH.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_char_p]
        self.rpa.performanceGetIdealIspH.restype = ctypes.c_double
        self.rpa.performancePrint.argtypes = [ctypes.c_void_p]

        self.rpa.chamberCreate.argtypes = [ctypes.c_void_p, ctypes.c_bool]
        self.rpa.chamberCreate.restype = ctypes.c_void_p
        self.rpa.chamberGetThrust.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        self.rpa.chamberGetThrust.restype = ctypes.c_double
        self.rpa.chamberGetMdot.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        self.rpa.chamberGetMdot.restype = ctypes.c_double
        self.rpa.chamberGetMdotOxidizer.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        self.rpa.chamberGetMdotOxidizer.restype = ctypes.c_double
        self.rpa.chamberGetMdotFuel.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        self.rpa.chamberGetMdotFuel.restype = ctypes.c_double
        self.rpa.chamberGetThroatDiameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        self.rpa.chamberGetThroatDiameter.restype = ctypes.c_double
        # self.rpa.nozzleCreate.argtypes = [ctypes.c_void_p, ctypes.c_bool]
        # self.rpa.nozzleCreate.restype = ctypes.c_void_p
        # self.rpa.nozzleGetExitDiameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        # self.rpa.nozzleGetExitDiameter.restype = ctypes.c_double

        self.rpa.nozzleStationCreate.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p, ctypes.c_int, ctypes.c_bool, ctypes.c_bool, ctypes.c_double, ctypes.c_char_p, ctypes.c_bool]
        self.rpa.nozzleStationCreate.restype = ctypes.c_void_p
        self.rpa.nozzleStationGetW.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        self.rpa.nozzleStationGetW.restype = ctypes.c_double

        resources = f'{path}/resources'
        self.rpa.initializeWithPath(resources.encode('utf-8'), None, 1)
        self.config = self.rpa.configFile()

    def sim_config(self, p_c, fuel, oxidizer, o_f_ratio, p_c_units='bar'):
        # fuel should be 'HTPB' and oxidizer should be 'N2O(L),298.15K'
        self.rpa.configFileCombustionChamberConditionsSetPressure(self.config, p_c, p_c_units.encode('utf-8'))
        self.rpa.configFilePropellantAddOxidizer(self.config, oxidizer.encode('utf-8'), 1.0, 0.0, 'bar'.encode('utf-8'), 0.0, 'K'.encode('utf-8'))
        self.rpa.configFilePropellantAddFuel(self.config, fuel.encode('utf-8'), 1, 0.0, 'bar'.encode('utf-8'), 0.0, 'K'.encode('utf-8'))
        self.rpa.configFilePropellantSetRatio(self.config, o_f_ratio, 'O/F'.encode('utf-8'))

    def model_chamber(self, thrust=None, m_dot=None, d_throat=None, l_star=None, contraction_angle=None, p_ambient=None):
        err = ValueError('One chamber condition is required!')
        if thrust is None and m_dot is None and d_throat is None:
            raise err
        if thrust is not None:
            if m_dot is not None or d_throat is not None:
                raise err
            if p_ambient is None:
                p_ambient = [1, 'bar']
            self.rpa.configFileEngineSizeSetThrust(self.config, thrust[0], thrust[1].encode('utf-8'))
            self.rpa.configFileEngineSizeSetAmbientPressure(self.config, p_ambient[0], p_ambient[1].encode('utf-8'))

        if m_dot is not None:
            if d_throat is not None:
                raise err
            self.rpa.configFileEngineSizeSetMdot(self.config, m_dot[0], m_dot[1].encode('utf-8'))

        if d_throat is not None:
            self.rpa.configFileEngineSizeSetThroatDiameter(self.config, d_throat[0], d_throat[1].encode('utf-8'))

        if l_star is not None:
            self.rpa.configFileEngineSizeSetChamberLength(self.config, l_star[0], l_star[1].encode('utf-8'), 1)

        if contraction_angle is not None:
            self.rpa.configFileEngineSizeSetContractionAngle(self.config, contraction_angle, 'degrees'.encode('utf-8'))

    def save_config_file(self, filename):
        self.rpa.configFileSaveAs(self.config, filename.encode('utf-8'))

    def model_nozzle_flow(self, exit_condition, exit_condition_value, contraction_ratio=None):
        '''
        Exit condition must be one of 'area ratio', 'pressure ratio', or 'exit pressure'
        'Pa'|'MPa'|'atm'|'bar'|'psi'|'psia'|'at'|'kg/cm2'
        '''
        self.rpa.configFileNozzleFlowOptionsSetCalculateNozzleFlow(self.config, True)
        if contraction_ratio is not None:
            self.rpa.configFileNozzleFlowOptionsSetNozzleInletConditions(self.config, 0, contraction_ratio)
        exit_conditions = {'area ratio': 0, 'pressure ratio': 1, 'exit pressure': 2}
        self.rpa.configFileNozzleFlowOptionsSetNozzleExitConditions(self.config, exit_conditions[exit_condition], exit_condition_value, None)

    def solve_performance(self):
        self.performance = self.rpa.performanceCreate(self.config, True, False)
        self.o_f_ratio = self.rpa.performanceGetOF(self.performance)
        self.isp = self.rpa.performanceGetDeliveredIsp(self.performance, 's'.encode('utf-8'), 0.0, 'bar'.encode('utf-8'), 0)
        self.throat = self.rpa.nozzleStationCreate(self.performance, 0.0, None, 2, False, False, 1, 'bar'.encode('utf-8'), False)
        self.throat_velocity = self.rpa.nozzleStationGetW(self.throat, 'm'.encode('utf-8'))

    def build_chamber(self,  thrust_units='N', m_dot_units='kg/s', throat_units='mm'):
        self.chamber = self.rpa.chamberCreate(self.performance, True)
        self.thrust = self.rpa.chamberGetThrust(self.chamber, thrust_units.encode('utf-8'))
        self.m_dot = self.rpa.chamberGetMdot(self.chamber, m_dot_units.encode('utf-8'))
        self.m_dot_ox = self.rpa.chamberGetMdotOxidizer(self.chamber, m_dot_units.encode('utf-8'))
        self.m_dot_fuel = self.rpa.chamberGetMdotFuel(self.chamber, m_dot_units.encode('utf-8'))
        self.d_throat = self.rpa.chamberGetThroatDiameter(self.chamber, throat_units.encode('utf-8'))
        # self.rpa.nozzleCreate.argtypes = [ctypes.c_void_p, ctypes.c_bool]
        # self.rpa.nozzleCreate.restype = ctypes.c_void_p
        # self.rpa.nozzleGetExitDiameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        # self.rpa.nozzleGetExitDiameter.restype = ctypes.c_double

# class HM9(RPASim):
#     """docstring for HM9."""
#
#     def __init__(self, m_dot_fuel, m_dot_ox, path='/Users/arnie/Downloads/rpa_lib'):
#         super().__init__(path)
#         self.o_f_ratio = m_dot_ox / m_dot_fuel
#         self.sim_config(p_c=40.0, fuel='HTPB', oxidizer='N2O(L),298.15K', o_f_ratio=self.o_f_ratio)
#         self.rebuild(m_dot_fuel, m_dot_ox)
#
    def rebuild(self, p_c, m_dot_fuel, m_dot_ox):
        self.o_f_ratio = m_dot_ox / m_dot_fuel
        self.rpa.configFileCombustionChamberConditionsSetPressure(self.config, p_c, 'bar'.encode('utf-8'))
        self.rpa.configFilePropellantSetRatio(self.config, self.o_f_ratio, 'O/F'.encode('utf-8'))
        self.solve_performance()
        self.build_chamber(throat_units='m')
