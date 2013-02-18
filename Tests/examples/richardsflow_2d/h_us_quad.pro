<?xml version="1.0"?>
<ogs6>
<coupling>
	<P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
		<out>PRESSURE</out>
		<out>SATURATION1</out>
		<problems>
			<M name="RICHARDS_FLOW" type="RICHARDS_FLOW">
				<out>PRESSURE1</out>
				<out>SATURATION1</out>
			</M>
		</problems>
	</P>
</coupling>
</ogs6>

