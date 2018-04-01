"""
RGBD Sensor factory
Author: Jeff Mahler
"""
from . import Kinect2Sensor

class RgbdSensorFactory:
    """ Factory class for Rgbd camera sensors. """

    @staticmethod
    def sensor(sensor_type, cfg):
        """ Creates a camera sensor of the specified type.

        Parameters
        ----------
        sensor_type : :obj:`str`
            the type of the sensor (real or virtual)
        cfg : :obj:`YamlConfig`
            dictionary of parameters for sensor initialization
        """
        sensor_type = sensor_type.lower()
        if sensor_type == 'kinect2':
            s = Kinect2Sensor(packet_pipeline_mode=cfg['pipeline_mode'],
                              device_num=cfg['device_num'],
                              frame=cfg['frame'])
        elif sensor_type == 'virtual_kinect2':
            pass
        else:
            raise ValueError('RGBD sensor type %s not supported' %(sensor_type))
        return s
