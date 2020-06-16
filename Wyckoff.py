import numpy as np
import re

class Position(object):
    def __init__(self,coordinate_functions):
        
        self.coordinate_functions = coordinate_functions

    @staticmethod
    def generate_from_string(string,base_string = '(0,0,0)'):
        
        positions_string = re.findall(r'\([^\)]+\)',string)
        base_positions_string = re.findall(r'\([^\)]+\)',base_string)

        coordinate_functions = list()
        for base_posistion_string in base_positions_string:
            base_coord_function = CoordinateFunction.generate_from_string(base_posistion_string)     
            for posistion_string in positions_string:
                coordinate_function = CoordinateFunction.generate_from_string(posistion_string) + base_coord_function
                coordinate_functions.append(coordinate_function)
        
        return Position(coordinate_functions)

    def __call__(self,**kwagrs):

        variable_values = np.zeros((3,))
        for var, value in kwagrs.items():
            variable_values['xyz'.index(var)] = value

        coordinates = np.zeros((len(self.coordinate_functions),3))
        for i, coordinate_function in enumerate(self.coordinate_functions):
            coordinate = coordinate_function(variable_values)
            coordinates[i,:] = coordinate

        return coordinates

class CoordinateFunction(object):
    
    def __init__(self,coord_matrix):
        
        self.coord_matrix = coord_matrix

    def __str__(self):

        return str(self.variables()) + '\n' + str(self.coord_matrix)

    def __add__(self,other):

        return self.__class__(self.coord_matrix + other.coord_matrix)

    def __call__(self,values):

        vector = np.ones((4,))
        vector[:3] = values

        return np.matmul(self.coord_matrix,vector)

    def variables(self):

        variables_string = np.array(['x','y','z'])
        variables = variables_string[np.any(self.coord_matrix[:,:3] != 0,axis = 0)]

        return variables

    @staticmethod
    def generate_from_string(string):
        coord_strings = re.split(r'[\(\),]',string)[1:4]
        coord_matrix = np.zeros((3,4))
        for i, coord_string in enumerate(coord_strings):
            expressions_string = re.findall(r'([-+]?)([xyz./0-9]+)',coord_string)
            for expression_string in expressions_string:
                if expression_string[0] is '-':
                    sign = -1
                else:
                    sign = +1

                if expression_string[1] in 'xyz':
                    j = 'xyz'.index(expression_string[1])
                    coord_matrix[i,j] = sign
                else:
                    coord_matrix[i,3] = eval(expression_string[1])

        return CoordinateFunction(coord_matrix)



        
