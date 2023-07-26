def update():
    """
    Download the newest version of files with solar indices. The files are updated daily.
    """
    # TODO: Save latest update date and implement autoupdating
    from urllib import request, error
    import os
    location = os.path.dirname(os.path.abspath(__file__))
    try:
        request.urlretrieve(
            'https://chain-new.chain-project.net/echaim_downloads/apf107.dat',
            os.path.join(location, f'data/index/apf107.dat')
        )
        request.urlretrieve(
            'https://chain-new.chain-project.net/echaim_downloads/ig_rz.dat',
            os.path.join(location, f'data/index/ig_rz.dat')
        )
        print('The index data was successfully updated!')
    except error.URLError:
        print('Something went wrong. Check internet connection.')
