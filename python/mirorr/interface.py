"""Mirorr symmetric rigid registration"""


import os
from nipype.interfaces.base import (
    CommandLine, CommandLineInputSpec, TraitedSpec, File,
    traits, isdefined)


class MirorrInputSpec(CommandLineInputSpec):
    moving_image = File(desc='Image to register', mandatory=True, exists=True,
                        argstr='--moving %s')
    moving_image_mask = File(
        desc='mask used to limit metric sampling region of the moving image',
        exists=True, argstr = '--moving-mask %s')
    fixed_image = File(desc='Image to register', mandatory=True, exists=True,
                       argstr='--fixed %s')
    fixed_image_mask = File(
        desc='mask used to limit metric sampling region of the fixed image',
        exists=True, argstr = '--fixed-mask %s')
    initial_moving_transform = File(desc='File with initial transformation',
                                    exists=True, argstr='--start-tfm %s')
    invert_initial_moving_transform = traits.Bool(
        desc='Whether to invert the transform.', argstr='--invert-start-tfm')
    transform = File(
        desc='Optimised transform (output name or already calculated)',
        argstr='--last-tfm %s', genfile=True)
    invert_output_transform = traits.Bool(
        desc='Whether to invert the transform.', argstr='--invert-last-tfm')
    transform_type = traits.Enum('rigid', 'affine', argstr='--tfm-type %s',
                                 usedefault=True)
    perform_registration = traits.Bool(
        True, desc='Perform the registration (or just resample, etc.).',
        argstr='--do-not-register')
    reoriented_fixed = File(
        desc='Filename to save fixed image in moving image space after '
        'registration. If rigid transform, only direction & origin updated, '
        'no resampling occurs',
        argstr='--save-reoriented-fixed %s', xor=['resampled_fixed'])
    reoriented_moving = File(
        desc='Filename to save moving image in fixed image space after '
        'registration. If rigid transform, only direction & origin updated, '
        'no resampling occurs',
        argstr='--save-reoriented-moving %s', xor=['resampled_moving'])
    resampled_fixed = File(
        desc='Filename to save resampled fixed image in moving image space '
        'after registration.',
        argstr='--save-resampled-fixed %s', xor=['reoriented_fixed'])
    resampled_moving = File(
        desc='Filename to save resampled moving image in fixed image space '
        'after registration.',
        argstr='--save-resampled-moving %s', xor=['reoriented_moving'])
    # Use ANTs-style eumerables for more homogeneous interface
    interpolation = traits.Enum(
        'Linear', 'NearestNeighbor', 'Sinc', 'BSpline',
        argstr='--final-interpolater %s')

    pyramid_start = traits.Int(
        1,
        desc='First pyramid level to be processed (-ve numbers count back '
        'from max).',
        argstr='--pyr-start %d')
    pyramid_end = traits.Int(
        desc='Last pyramid level to be processed (-ve numbers count back '
        'from max).',
        argstr='--pyr-end %d')
    number_of_pyramid_levels = traits.Int(
        1024,
        desc='Number of levels in the pyramid (default is to create as many '
        'levels as maxDataSize > 32)',
        argstr='--pyr-num')
    number_of_iterations = traits.Int(
        desc='Number of iterations. Set negative to halve iterations each '
        'pyramid level.',
        argstr='--iterations %d')
    metric = traits.Enum(
        'NC', 'SD', 'CR', 'MI',
        desc='metric used: normalized correlation (NC), sum of squared '
        'differences (SD), correlation ratio (CR), non-parametric window '
        'mutual information (MI).',
        argstr='--blockmetric %s')
    block_width = traits.Int(
        desc='Block size in pixels', argstr='--blockwidth %d')

class MirorrOutputSpec(TraitedSpec):
    forward_transforms = traits.List(
        File(exists=True),
        desc='List of output transforms for forward registration')
    reverse_transforms = traits.List(
        File(exists=True),
        desc='List of output transforms for reverse registration')
    forward_invert_flags = traits.List(
        traits.Bool(),
        desc='List of flags corresponding to the forward transforms')
    reverse_invert_flags = traits.List(
        traits.Bool(),
        desc='List of flags corresponding to the reverse transforms')
    warped_image = File(desc='Alias for resampled_moving')
    inverse_warped_image = File(desc='Alias for resampled_fixed')
    reoriented_fixed = File(desc='Output reoriented fixed image.')
    reoriented_moving = File(desc='Output reoriented moving image.')
    resampled_fixed = File(desc='Output resampled fixed image.')
    resampled_moving = File(desc='Output resampled moving image.')

class Mirorr(CommandLine):
    """Mirror symmetric rigid registration."""
    _cmd = 'mirorr'
    input_spec = MirorrInputSpec
    output_spec = MirorrOutputSpec

    def _list_outputs(self):
        outputs = self._outputs().get()

        # get output transforms
        outputs['reverse_transforms'] = outputs['forward_transforms'] = \
            [self._gen_filename('transform')]
        outputs['forward_invert_flags'] = [False]
        outputs['reverse_invert_flags'] = [True]

        # get output images if they were calculated
        for output_name in ('reoriented_fixed', 'reoriented_moving',
                            'resampled_fixed', 'resampled_moving'):
            input_val = getattr(self.inputs, output_name)
            if isdefined(input_val):
                outputs[output_name] = os.path.abspath(input_val)

        # aliases
        if 'resampled_fixed' in outputs:
            outputs['warped_image'] = outputs['resampled_moving']
        if 'resampled_moving' in outputs:
            outputs['inverse_warped_image'] = outputs['resampled_fixed']

        return outputs

    def _gen_filename(self, name):
        if name == 'transform':
            if isdefined(self.inputs.transform):
                return self.inputs.transform
            else:
                return os.path.abspath(
                    '{}.tfm'.format(self.inputs.transform_type))
        else:
            return super()._gen_filename(name)

    def _format_arg(self, name, trait_spec, value):
        """milxOrientFilterUsingDirectionCosines expects `orientation` to be a
        number, the index in `orientations`."""
        if name == 'perform_registration':
            # invert to use --do-not-register flag if not perform_registration
            value = not value
        if name == 'interpolation':
            # Go from ANTs-style definition to Mirorr's
            value = dict(
                Linear='linear',
                NearestNeighbor='nn',
                BSpline='bspline',
                Sinc='sinc')[value]
        if name == 'metric':
            # lower case metric names
            value = value.lower()
        return super()._format_arg(name, trait_spec, value)
