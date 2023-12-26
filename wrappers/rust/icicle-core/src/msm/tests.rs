use super::MSM;
use crate::curve::{CurveConfig, Projective};
use crate::field::FieldConfig;
use crate::traits::GenerateRandom;
use icicle_cuda_runtime::{memory::DeviceSlice, stream::CudaStream};

#[cfg(feature = "arkworks")]
use crate::traits::ArkConvertible;
#[cfg(feature = "arkworks")]
use ark_ec::models::CurveConfig as ArkCurveConfig;
#[cfg(feature = "arkworks")]
use ark_ec::VariableBaseMSM;

pub fn check_msm<C: CurveConfig + MSM<C>, ScalarConfig: FieldConfig>()
where
    ScalarConfig: GenerateRandom<C::ScalarField>,
    C::ScalarField: ArkConvertible,
    C::ScalarField: ArkConvertible<ArkEquivalent = <C::ArkSWConfig as ArkCurveConfig>::ScalarField>,
    C::BaseField: ArkConvertible,
    C::BaseField: ArkConvertible<ArkEquivalent = <C::ArkSWConfig as ArkCurveConfig>::BaseField>,
{
    let test_sizes = [1000, 1 << 18];
    let mut msm_results = DeviceSlice::cuda_malloc(1).unwrap();
    for test_size in test_sizes {
        let points = C::generate_random_affine_points(test_size);
        let scalars = ScalarConfig::generate_random(test_size);
        let stream = CudaStream::create().unwrap();

        let mut cfg = C::get_default_msm_config();
        cfg.ctx
            .stream = &stream;
        cfg.is_async = true;
        cfg.are_results_on_device = true;
        C::msm(&scalars, &points, cfg, &mut msm_results.as_slice_mut()).unwrap();

        let mut msm_host_result = vec![Projective::<C>::zero(); 1];
        msm_results
            .copy_to_host(&mut msm_host_result[..])
            .unwrap();
        stream
            .synchronize()
            .unwrap();
        stream
            .destroy()
            .unwrap();

        let point_r_ark: Vec<_> = points
            .iter()
            .map(|x| x.to_ark())
            .collect();
        let scalars_r_ark: Vec<_> = scalars
            .iter()
            .map(|x| x.to_ark())
            .collect();
        let msm_result_ark: ark_ec::models::short_weierstrass::Projective<C::ArkSWConfig> =
            VariableBaseMSM::msm(&point_r_ark, &scalars_r_ark).unwrap();
        assert_eq!(msm_host_result[0].to_ark(), msm_result_ark);
    }
}
