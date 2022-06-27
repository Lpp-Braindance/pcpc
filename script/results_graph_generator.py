from matplotlib import pyplot
from matplotlib.axes import Axes

# x-axis values 
x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10 , 11 , 12] 
# y-axis values 

pyplot.title("Strong scalability")
y = [ 173.624 ,123.140 ,84.922 ,61.565 ,53.424 ,48.457 ,43.463 ,37.239 ,35.890 ,30.699 ,22.441 ,20.770 ] 

# pyplot.title("Weak scalability")
# y = [ 1.203 ,3.434 ,5.130 ,6.888 ,8.577 ,10.413 ,12.020 ,13.807 ,15.519 ,17.228 ,18.849 ,20.592 ] 


pyplot.xlabel("processes")
pyplot.ylabel("Time[sec]")

pyplot.xlim([-1, 13])

pyplot.plot(x, y)
pyplot.grid()
pyplot.show() 


