<script>
  import svelteLogo from './assets/svelte.svg'
  import Counter       from './lib/Counter.svelte'
  import FritexpToyOne from './lib/FritexpToyOne.svelte'
  import fritexp       from './lib/fritexp_szekeres.js'
  
  let xvalue; 
  let uvalue;
  let answer;
  $: answer = fritexp(uvalue, xvalue);
</script>

<main>
 <h1>Fractionally Iterated Exponential</h1>
 <FritexpToyOne />
 <p>exp_1/2(0) = { fritexp(0.5, 0) }</p>
 <label for="xvalue">x</label>
 <input type="number" id="xvalue" name="x" bind:value={xvalue}>
 <label for="uvalue">u</label>
 <input type="number" id="uvalue" name="u" bind:value={uvalue}>
 <p>exp_({uvalue})({xvalue}) = {answer}</p>
</main>

<style>
  
</style>
