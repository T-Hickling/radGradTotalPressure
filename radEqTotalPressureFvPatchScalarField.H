/*---------------------------------------------------------------------------*\
Class
    Foam::radEqTotalPressureFvPatchScalarField

Description
    This boundary condition provides a total pressure condition for patches
    in a rotating frame.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | none
        psi          | compressibility field name | no       | none
        gamma        | ratio of specific heats (Cp/Cv) | yes |
        p0           | static pressure reference | yes       |
        omega        | angular velocity of the frame [rad/s] | yes    | none
        dp0dr        | radial total pressure gradient | yes | none
        rRef             | reference radial location for p0   | none
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            rotEqTotalPressure;
        U               U;
        phi             phi;
        rho             rho;
        psi             psi;
        gamma           1.4;
        p0              uniform 1e5;
        omega           (0 0 100);
    }
    \endverbatim

    The \c omega entry is a Function1 type, able to describe time varying
    functions.

See also
    Foam::totalPressureFvPatchScalarField
    Foam::Function1Types

SourceFiles
    radEqTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rotEqTotalPressureFvPatchScalarField_H
#define rotEqTotalPressureFvPatchScalarField_H

#include "totalPressureFvPatchScalarField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class radEqTotalPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/


class radEqTotalPressureFvPatchScalarField
:
    public totalPressureFvPatchScalarField
{
    // Private Data

        //- Angular velocity of the frame
        autoPtr<Function1<vector>> omega_;

        //- Radial gradient of p0
        autoPtr<Function1<scalar>> dp0dr_;

        //- Radial location of p0
        scalar rRef_;


public:

    //- Runtime type information
    TypeName("radEqTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        radEqTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        radEqTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given radEqTotalPressureFvPatchScalarField
        //  onto a new patch
        radEqTotalPressureFvPatchScalarField
        (
            const radEqTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        radEqTotalPressureFvPatchScalarField
        (
            const radEqTotalPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        tmp<fvPatchScalarField> clone() const override
        {
            return tmp<fvPatchScalarField>
            (
                new radEqTotalPressureFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        radEqTotalPressureFvPatchScalarField
        (
            const radEqTotalPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchScalarField>
            (
                new radEqTotalPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            void updateCoeffs() override ;


        //- Write
        void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
