

from postqe import get_charge

charge = get_charge(prefix='Si', schema='../../schemas/qes.xsd')

#figure = charge.plot(x0=(0, 0, 0), e1=(1, 0, 0), nx=10, dim=1, ifmagn='total',
#             plot_file='postqeiflag0out0', method='spherical', format='gnuplot')
#figure = charge.plot(x0=(0, 0, 0), e1=(1, 0, 0), nx=10, dim=1, ifmagn='total',
#             plot_file='postqeiflag1out0', method='FFT', format='gnuplot')
#figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, dim=2, ifmagn='total',
#             plot_file='postqeiflag2out7', method='FFT', format='gnuplot')
#figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, dim=2, ifmagn='total',
#             plot_file='postqeiflag2out1', method='FFT', format='contour.x')
#figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, dim=2, ifmagn='total',
#             plot_file='postqeiflag2out2', method='FFT', format='plotrho.x')
#figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, dim=2, ifmagn='total',
#             plot_file='postqeiflag2out3', method='FFT', format='xsf')
#figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, radius = 10, dim=2, ifmagn='total',
#             plot_file='postqeiflag4out7', method='polar', format='gnuplot')
#figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, radius = 10, dim=2, ifmagn='total',
#             plot_file='postqeiflag4out1', method='polar', format='contour.x')
#figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, e3 = (0., 0., 1.), nz=10,
#                dim=3, ifmagn='total', plot_file='postqeiflag3out3', method='FFT', format='xsf')
figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, e3 = (0., 0., 1.), nz=10,
                dim=3, ifmagn='total', plot_file='postqeiflag3out6', method='FFT', format='cube')

figure.savefig("figure_1.png")
figure.savefig("figure_1.pdf", format='pdf')