<script>
  import svelteLogo from './assets/svelte.svg'
  import Counter       from './lib/Counter.svelte'
  import FritexpToyOne from './lib/FritexpToyOne.svelte'
  import {fritexp, fritexp_function_factory}  from './lib/fritexp_szekeres.js'
  
  let xvalue; 
  let uvalue;
  let answer;
  $: answer = fritexp(uvalue, xvalue);
  
  let calculator;
  let calculator_answer;
  $: calculator = fritexp_function_factory(uvalue);
  $: calculator_answer = calculator(xvalue);
</script>

<main>
 <h1>Fractionally Iterated Exponential</h1>
 <FritexpToyOne />
 <p>exp_1/2(0) = { fritexp(0.5, 0) }</p>
 <label for="xvalue">x=</label>
 <input type="number" id="xvalue" name="x" bind:value={xvalue}>
 <span style="width:13rem;">, </span>
 <label for="uvalue">u=</label>
 <input type="number" id="uvalue" name="u" bind:value={uvalue}>
 <p>exp_({uvalue})({xvalue}) = {answer},  calculator: {calculator_answer}</p>
</main>

<style>
  
</style>
