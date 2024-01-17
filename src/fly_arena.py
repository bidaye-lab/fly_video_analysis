from scipy.io import loadmat
import h5py


def load_track(matlab_file):
    # https://github.com/jstaf/fly_tracker


        # load tracking data
    try:
        # matlab file pre 7.3
        m = loadmat(matlab_file, squeeze_me=True, struct_as_record=False)
        data = vars(m['trk'])['data'][:, :, [0, 1]]
    except NotImplementedError:
        # matlab file since 7.3
        with h5py.File(matlab_file, 'r') as f:
            data = f['trk']['data'][()].T[:, :, [0, 1]]

    return data
    