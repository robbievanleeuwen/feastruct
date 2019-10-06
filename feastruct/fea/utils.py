import numpy as np


def gauss_points(el_type, n):
    """Returns the Gaussian weights and locations for *n* point Gaussian integration of a finite
    element. Refer to xxx for a list of the element types.

    :param string el_type: String describing the element type
    :param int n: Number of Gauss points
    :returns: The integration weights *(n x 1)* and an *(n x i)* matrix consisting of the values of
        the *i* shape functions for *n* Gauss points
    :rtype: tuple(list[float], :class:`numpy.ndarray`)
    """

    if el_type == 'Tri6':
        # one point gaussian integration
        if n == 1:
            weights = [1]
            gps = np.array([[1.0 / 3, 1.0 / 3, 1.0 / 3]])

        # three point gaussian integration
        elif n == 3:
            weights = [1.0 / 3, 1.0 / 3, 1.0 / 3]
            gps = np.array([
                [2.0 / 3, 1.0 / 6, 1.0 / 6],
                [1.0 / 6, 2.0 / 3, 1.0 / 6],
                [1.0 / 6, 1.0 / 6, 2.0 / 3]
            ])

        # six point gaussian integration
        elif n == 6:
            g1 = 1.0 / 18 * (8 - np.sqrt(10) + np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
            g2 = 1.0 / 18 * (8 - np.sqrt(10) - np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
            w1 = (620 + np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
            w2 = (620 - np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720

            weights = [w2, w2, w2, w1, w1, w1]
            gps = np.array([
                [1 - 2 * g2, g2, g2],
                [g2, 1 - 2 * g2, g2],
                [g2, g2, 1 - 2 * g2],
                [g1, g1, 1 - 2 * g1],
                [1 - 2 * g1, g1, g1],
                [g1, 1 - 2 * g1, g1]
            ])

    return (weights, gps)


def shape_function(el_type, coords, gp):
    """Computes shape functions, shape function derivatives and the determinant of the Jacobian
    matrix for a number of different finite elements at a given Gauss point. Refer to xxx for a
    list of the element types.

    :param string el_type: String describing the element type
    :param coords: Global coordinates of the element nodes *(n x 3)*, where *n* is the number of
        nodes
    :type coords: :class:`numpy.ndarray`
    :param gp: Isoparametric location of the Gauss point
    :type gp: :class:`numpy.ndarray`
    :returns: The value of the shape functions *N(i)* at the given Gauss point *(1 x n)*, the
        derivative of the shape functions in the j-th global direction *B(i,j)* *(3 x n)* and the
        determinant of the Jacobian matrix *j*
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`, float)
    """

    if el_type == 'Tri6':
        # location of isoparametric co-ordinates for each Gauss point
        eta = gp[0]
        xi = gp[1]
        zeta = gp[2]

        # value of the shape functions
        N = np.array([
            eta * (2 * eta - 1),
            xi * (2 * xi - 1),
            zeta * (2 * zeta - 1),
            4 * eta * xi,
            4 * xi * zeta,
            4 * eta * zeta
        ])

        # derivatives of the sf wrt the isoparametric co-ordinates
        B_iso = np.array([
            [4 * eta - 1, 0, 0, 4 * xi, 0, 4 * zeta],
            [0, 4 * xi - 1, 0, 4 * eta, 4 * zeta, 0],
            [0, 0, 4 * zeta - 1, 0, 4 * xi, 4 * eta]
        ])

        # form Jacobian matrix
        J_upper = np.array([[1, 1, 1]])
        J_lower = np.dot(coords, np.transpose(B_iso))
        J = np.vstack((J_upper, J_lower))

        # calculate the jacobian
        j = 0.5 * np.linalg.det(J)

        # cacluate the P matrix
        P = np.dot(np.linalg.inv(J), np.array([[0, 0], [1, 0], [0, 1]]))

        # calculate the B matrix in terms of cartesian co-ordinates
        B = np.transpose(np.dot(np.transpose(B_iso), P))

    return (N, B, j)
