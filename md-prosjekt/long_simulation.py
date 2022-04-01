from click import style
from box import *
from n_atom_sim import *

figpath = os.path.join(os.getcwd(),'figures')

md = MD(r0=lattice(n=10),
        n=4000,
        bound=True,
        rc=3,
        test=True)
md.set_inital_velocities(T=180)
t, x, v = md.solve(10, 0.01)

A = md.vac()
plt.plot(t, A)
plt.axhline(0, linestyle='--', alpha=0.8, color='k')
plt.xlabel('t*')
plt.ylabel('Velocity autocorrelation')
plt.grid()
plt.savefig(os.path.join(os.getcwd(),'figures/long_sim_vac.pdf'))
plt.clf()
# plt.show()

D = md.diffusion_coefficient()
plt.plot(t, D)
plt.xlabel('t*')
plt.ylabel('Diffusion coefficient')
plt.grid()
plt.savefig(os.path.join(os.getcwd(),'figures/long_sim_diff.pdf'))
plt.clf()
# plt.show()

msd, t = md.msd()
plt.plot(t, msd)
plt.xlabel('t*')
plt.ylabel('Mean square displacement')
plt.grid()
plt.savefig(os.path.join(os.getcwd(),'figures/long_sim_msd.pdf'))
plt.clf()
# plt.show()

rdf, bins = md.rdf(200)
plt.plot(bins, rdf)
plt.xlabel('r*')
plt.ylabel('Radial distribution')
plt.axhline(1, linestyle='--', alpha=0.8, color='k')
plt.savefig(os.path.join(os.getcwd(),'figures/long_sim_rdf.pdf'))
plt.clf()
# plt.show()
