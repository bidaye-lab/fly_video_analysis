import numpy as np


def adjust_crop(vid, x, y, d):
    """Check if crop area is within video and adjust accordingly


    Parameters
    ----------
    vid : numpy.array
        video
    x : int
        x coordinate of the center in the cropped video
    y : int
        x coordinate of the center in the cropped video
    d : int
        dimension in pixel of the cropped video

    Returns
    -------
    x : int
        adjusted x coordinate of the center in the cropped video
    y : int
        adjusted y coordinate of the center in the cropped video
    """

    _, res_y, res_x = vid.shape  # get video resolution

    # adjust lower bound
    x, y = np.max([d // 2, x]), np.max([d // 2, y])
    # adjust upper bound
    x, y = np.min([res_x - (d // 2), x]), np.min([res_y - (d // 2), y])

    return x, y
