#include "radEqTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radEqTotalPressureFvPatchScalarField::
radEqTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(p, iF),
    omega_(),
    dp0dr_(),
    rRef_(0.0)
{}


Foam::radEqTotalPressureFvPatchScalarField::
radEqTotalPressureFvPatchScalarField
(
    const radEqTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(ptf, p, iF, mapper),
    omega_(ptf.omega_, false),
    dp0dr_(ptf.dp0dr_, false),
    rRef_(ptf.rRef_)
{}


Foam::radEqTotalPressureFvPatchScalarField::
radEqTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF, dict),
    omega_(Function1<vector>::New("omega", dict)),
    dp0dr_(Function1<scalar>::New("dp0dr", dict)),
    rRef_(readScalar(dict.lookup("rRef")))
{}


Foam::radEqTotalPressureFvPatchScalarField::
radEqTotalPressureFvPatchScalarField
(
    const radEqTotalPressureFvPatchScalarField& retppsf
)
:
    totalPressureFvPatchScalarField(retppsf),
    omega_(retppsf.omega_, false),
    dp0dr_(retppsf.dp0dr_, false),
    rRef_(retppsf.rRef_)
{}


Foam::radEqTotalPressureFvPatchScalarField::
radEqTotalPressureFvPatchScalarField
(
    const radEqTotalPressureFvPatchScalarField& retppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(retppsf, iF),
    omega_(retppsf.omega_, false),
    dp0dr_(retppsf.dp0dr_, false),
    rRef_(retppsf.rRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radEqTotalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    const vector om = omega_->value(t);
    const scalar r0 = rRef_;
    const scalar pGrad = dp0dr_->value(t);

    vector axisHat = om/mag(om);
    tmp<vectorField> rotationVelocity =
        om ^ (patch().Cf() - axisHat*(axisHat & patch().Cf()));

    tmp<scalarField> r = mag(patch().Cf() ^ axisHat);

    const scalarField pCorr
    (
        pGrad*(r - r0)
    );

    const vectorField Up
    (
        patch().lookupPatchField<volVectorField, vector>(UName())
      + rotationVelocity
    );

    totalPressureFvPatchScalarField::updateCoeffs(p0() + pCorr, Up);
}


void Foam::radEqTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);
    writeEntry(os, omega_());
    writeEntry(os, dp0dr_());
    writeEntry(os, "rRef", rRef_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
    fvPatchScalarField,
    radEqTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //