<script>
import  ValueSlider   from './ValueSlider.svelte'
let order_iter = 0.0;
let order_minval = -5.2;
let order_maxval =  5.2;

</script>



<ValueSlider label="Order of Iteration" minval={order_minval} maxval={order_maxval} bind:value={order_iter} />
<FunctionPlot />

<style>
</style>
