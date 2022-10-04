def update():
    from urllib import request, error
    import os
    location = os.path.dirname(os.path.abspath(__file__))
    try:
        for ver in ['16', '20']:
            request.urlretrieve(
                'https://chain-new.chain-project.net/echaim_downloads/apf107.dat',
                os.path.join(location, f'data/data{ver}/index/apf107.dat')
            )
            request.urlretrieve(
                'https://chain-new.chain-project.net/echaim_downloads/ig_rz.dat',
                os.path.join(location, f'data/data{ver}/index/ig_rz.dat')
            )
        print('The index data was successfully updated!')
    except error.URLError:
        print('Something went wrong. Check internet connection.')
