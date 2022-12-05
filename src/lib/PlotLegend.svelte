<script>
export let items = []

function addItem(name)  {
    let newcolor = "#22aa00";
    items.append( {color: newcolor, name: name} )
    return newcolor; 
}



function removeAllItems()  {
}
    
$: console.log(items)   // redraw legend
</script>



<div class="border">
<p>Pretend Legend</p>
<ul>
<li>----- blah</li>
<li>===== xyz</li>
<li>..... zzzz</li>
<li>ooooo more</li>
</ul>
</div>


<style>
.border { border: solid 4px #eee; border-radius:9px; }

</style>
