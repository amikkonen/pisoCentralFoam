fvVectorMatrix UrelEqn
(
    fvm::ddt(rho,Urel)
    +
    fvm::div(phiPos,Urel)
    +
    fvm::div(phiNeg,Urel)
    +
    rho*SRF->Su()
    +
    turbulence->divDevRhoReff(Urel)
);

if (pimple.momentumPredictor())
{
    solve
    (
	UrelEqn == -gradp
    );
}
