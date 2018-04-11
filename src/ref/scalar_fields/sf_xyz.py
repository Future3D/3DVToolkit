import numpy as np
from ..geometry.coord_systems import (
    cartesian_to_spherical,
    cartesian_to_cylindrical
)
from ..ransac import single_fit, RANSAC_MODELS, RANSAC_SAMPLERS
from .base import ScalarField


class ScalarField_XYZ(ScalarField):
    def extract_info(self):
        self.points = self.pyntcloud.xyz


class PlaneFit(ScalarField_XYZ):
    """
    Get which points belong to the best RansacPlane found.
    """

    def __init__(self, pyntcloud, max_dist=1e-4, max_iterations=100, n_inliers_to_stop=None):
        self.model = RANSAC_MODELS["plane"]
        self.sampler = RANSAC_SAMPLERS["random"]
        self.name = "is_plane"
        self.model_kwargs = {"max_dist": max_dist}
        self.max_iterations = max_iterations
        self.n_inliers_to_stop = n_inliers_to_stop

        super().__init__(pyntcloud)

    def compute(self):
        inliers = single_fit(self.points, self.model, self.sampler,
                             model_kwargs=self.model_kwargs,
                             max_iterations=self.max_iterations,
                             n_inliers_to_stop=self.n_inliers_to_stop)
        self.to_be_added[self.name] = inliers.astype(np.uint8)


class SphereFit(ScalarField_XYZ):
    """
    Get which points belong to the best RansacSphere found.
    """

    def __init__(self, pyntcloud, max_dist=1e-4, max_iterations=100, n_inliers_to_stop=None):
        self.model = RANSAC_MODELS["sphere"]
        self.sampler = RANSAC_SAMPLERS["random"]
        self.name = "is_sphere"
        self.model_kwargs = {"max_dist": max_dist}
        self.max_iterations = max_iterations
        self.n_inliers_to_stop = n_inliers_to_stop

        super().__init__(pyntcloud)

    def compute(self):
        inliers = single_fit(self.points, self.model, self.sampler,
                             model_kwargs=self.model_kwargs,
                             max_iterations=self.max_iterations,
                             n_inliers_to_stop=self.n_inliers_to_stop)
        self.to_be_added[self.name] = inliers.astype(np.uint8)


class CustomFit(ScalarField_XYZ):
    """
    Fit using custom model and sampler.
    """

    def __init__(self, pyntcloud, model, sampler, name, model_kwargs={},
                 sampler_kwargs={}, max_iterations=100, n_inliers_to_stop=None):
        self.model = model
        self.sampler = sampler
        self.name = name
        self.model_kwargs = model_kwargs
        self.sampler_kwargs = sampler_kwargs
        self.max_iterations = max_iterations
        super().__init__(pyntcloud)

    def compute(self):
        inliers = single_fit(self.points, self.model, self.sampler,
                             model_kwargs=self.model_kwargs,
                             max_iterations=self.max_iterations,
                             n_inliers_to_stop=self.n_inliers_to_stop)
        self.to_be_added[self.name] = inliers.astype(np.uint8)


class SphericalCoordinates(ScalarField_XYZ):
    """
    Get radial, azimuthal and polar values from x, y, z coordinates.
    """

    def __init__(self, pyntcloud, degrees=True):
        self.degrees = degrees
        super().__init__(pyntcloud)

    def compute(self):
        radial, theta, phi = cartesian_to_spherical(
            self.points, degrees=self.degrees)

        self.to_be_added["radial"] = radial
        self.to_be_added["polar"] = theta
        self.to_be_added["azimuthal"] = phi


class CylindricalCoordinates(ScalarField_XYZ):
    """
    Get ro and phi values from x, y, z coordinates.
    The z value in cylindrical coordinates remain unchanged.
    """

    def __init__(self, pyntcloud, degrees=True):
        self.degrees = degrees
        super().__init__(pyntcloud)

    def compute(self):
        ro, phi, z = cartesian_to_cylindrical(
            self.points, degrees=self.degrees)

        self.to_be_added["ro"] = ro
        self.to_be_added["phi"] = phi