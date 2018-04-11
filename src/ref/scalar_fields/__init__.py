from .sf_eigenvalues import (
    Anisotropy,
    Curvature,
    Eigenentropy,
    EigenSum,
    Linearity,
    Omnivariance,
    Planarity,
    Sphericity
)
from .sf_kneighbors import (
    EigenDecomposition,
    EigenValues,
    Normals,
)
from .sf_normals import (
    InclinationDegrees,
    InclinationRadians,
    OrientationDegrees,
    OrientationRadians
)
from .sf_rgb import (
    HueSaturationValue,
    RelativeLuminance,
    RGBIntensity
)
from .sf_voxelgrid import (
    VoxelN,
    VoxelX,
    VoxelY,
    VoxelZ,
    EuclideanClusters
)
from .sf_xyz import (
    PlaneFit,
    SphereFit,
    CustomFit,
    SphericalCoordinates,
    CylindricalCoordinates
)
ALL_SF = {
    # Eigenvalues
    'anisotropy': Anisotropy,
    'curvature': Curvature,
    'eigenentropy': Eigenentropy,
    'eigen_sum': EigenSum,
    'linearity': Linearity,
    'omnivariance': Omnivariance,
    'planarity': Planarity,
    'sphericity': Sphericity,
    # Kneighbors
    'eigen_decomposition': EigenDecomposition,
    'eigen_values': EigenValues,
    'normals': Normals,
    # Normals
    'inclination_deg': InclinationDegrees,
    'inclination_rad': InclinationRadians,
    'orientation_deg': OrientationDegrees,
    'orientation_rad': OrientationRadians,
    # RGB
    'hsv': HueSaturationValue,
    'relative_luminance': RelativeLuminance,
    'rgb_intensity': RGBIntensity,
    # Voxelgrid
    'voxel_n': VoxelN,
    'voxel_x': VoxelX,
    'voxel_y': VoxelY,
    'voxel_z': VoxelZ,
    'euclidean_clusters': EuclideanClusters,
    # XYZ
    'custom_fit': CustomFit,
    'plane_fit': PlaneFit,
    'sphere_fit': SphereFit,
    'spherical_coords': SphericalCoordinates,
    'cylindrical_coords': CylindricalCoordinates
}
