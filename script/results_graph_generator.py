from matplotlib import pyplot
from matplotlib.axes import Axes

# x-axis values 
x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10 , 11 , 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24] 
# y-axis values 

pyplot.title("Strong scalability")
y = [ 140.244 , 78.054, 88.75 ] 

# pyplot.title("Weak scalability")
# y = [ 1.203 ,3.434 ,5.130 ,6.888 ,8.577 ,10.413 ,12.020 ,13.807 ,15.519 ,17.228 ,18.849 ,20.592 ] 


pyplot.xlabel("processes")
pyplot.ylabel("Time[sec]")

pyplot.xlim([-1, 13])

pyplot.plot(x, y)
pyplot.grid()
pyplot.show() 


