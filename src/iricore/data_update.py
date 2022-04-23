def update():
    from urllib import request, error
    import os
    location = os.path.dirname(os.path.abspath(__file__))
    try:
        request.urlretrieve('http://irimodel.org/indices/apf107.dat', os.path.join(location, 'data/index/apf107.dat'))
        request.urlretrieve('http://irimodel.org/indices/ig_rz.dat', os.path.join(location, 'data/index/ig_rz.dat'))
        print('The index data was successfully updated!')
    except error.URLError:
        print('Something went wrong. Check internet connection.')
