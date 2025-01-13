import click
from astropy import coordinates
from ne2001 import density, ne_io

@click.group('ne2001')
def main():
    pass

@main.command()
@click.argument('ra', type=float)
@click.argument('dec', type=float)
@click.option('--distance', type=float, default=10)
def mwdm(ra, dec, distance):
    co = coordinates.SkyCoord(ra, dec, unit='deg')
    # or pass sexagesimal
    #    co = coordinates.SkyCoord(rastr, decstr, unit=(units.hourangle, units.deg))

    ne = density.ElectronDensity(**ne_io.Params())
    dm = ne.DM(co.galactic.l, co.galactic.b, distance)
    
    print(f'For (RA, Dec) = ({ra}, {dec}), (l, b) = ({co.galactic.l}, {co.galactic.b}), DM={dm} pc/cm3')
