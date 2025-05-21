load("utils1.6.sage")

def rotation_matrix(axis, theta):
    """
    Return the 3x3 rotation matrix for a rotation of angle theta (in radians)
    about the given axis. The axis is specified as an iterable of 3 numbers.
    """
    # Convert the axis to a normalized vector
    axis = vector(axis).normalized()
    x, y, z = axis
    c = cos(theta)
    s = sin(theta)
    t = 1 - c  # common factor

    # Construct the rotation matrix using Rodrigues' formula
    R = matrix(
        [
            [c + x**2 * t, x * y * t - z * s, x * z * t + y * s],
            [y * x * t + z * s, c + y**2 * t, y * z * t - x * s],
            [z * x * t - y * s, z * y * t + x * s, c + z**2 * t],
        ]
    )
    return R
