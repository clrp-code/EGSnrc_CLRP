:start geometry definition:

    :start geometry:
        name = phantom
        library = egs_rz

        number of shells = 2 2 2
        shell thickness =  0.5 1 2

        first plane = -4
        number of slabs = 2 8 2
        slab thickness  = 1 0.5 1

        :start media input:
            media = M0
        :stop media input:
    :stop geometry:

    simulation geometry = phantom

:stop geometry definition:

