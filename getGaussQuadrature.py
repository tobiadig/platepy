import numpy as np
def getGaussQuadrature(shape, nPoints):
    gaussQuadrature={'linear': {1:{'points': np.array([[0,0]]),
                                    'weights': np.array([2])},

                                    2:{'points': np.array([[-1*np.sqrt(3),0],
                                                            [1*np.sqrt(3),0]]),
                                    'weights': np.array([1,1])}},
        
                    'triangular': {1:{'points': np.array([[1/3, 1/3]]),
                                    'weights': np.array([1/2])},

                                    3:{'points': np.array([[1/6, 1/6],
                                                            [2/3,   1/6],
                                                            [1/6, 2/3]]),
                                    'weights': np.array([1/6, 1/6, 1/6])},

                                    4:{'points': np.array([[1/3, 1/3],
                                                            [1/5, 1/5],
                                                            [3/5, 1/5],
                                                            [1/5, 3/5]]),
                                    'weights': np.array([-27/96, 25/96, 25/96, 25/96])}},

                    'rectangular':{1:{'points': np.array([[0, 0]]),
                                    'weights': np.array([4])},

                                    4:{'points': np.array([[-1/np.sqrt(3), -1/np.sqrt(3)],
                                                            [1/np.sqrt(3),   -1/np.sqrt(3)],
                                                            [1/np.sqrt(3),1/np.sqrt(3)],
                                                            [-1/np.sqrt(3), 1/np.sqrt(3)]]),
                                    'weights': np.array([1, 1, 1, 1])},

                                    9:{'points': np.array([[-np.sqrt(3/5), -np.sqrt(3/5)],
                                                            [ 0      , -np.sqrt(3/5)],
                                                            [+np.sqrt(3/5), -np.sqrt(3/5)],
                                                            [-np.sqrt(3/5),    0    ],
                                                            [   0    ,    0    ] ,
                                                            [+np.sqrt(3/5),    0    ],
                                                            [-np.sqrt(3/5), +np.sqrt(3/5)],
                                                            [   0    , +np.sqrt(3/5)],
                                                            [+np.sqrt(3/5), +np.sqrt(3/5)]]),
                                    'weights': np.array([25, 40, 25, 40, 64, 40, 25, 40, 25])/81}}}
    return gaussQuadrature[shape][nPoints]['points'], gaussQuadrature[shape][nPoints]['weights']