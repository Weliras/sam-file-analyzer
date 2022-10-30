function RandomColor () {
    var r = Math.floor(Math.random() * 255);
    var g = Math.floor(Math.random() * 255);
    var b = Math.floor(Math.random() * 255);
    return "rgb(" + r + "," + g + "," + b + ")";
}


const virus_coverage_json = JSON.parse(virus_coverage)
const genes_coverage_json = JSON.parse(genes_coverage)


// Virus coverage graph

let virus_coverage_labels = []
let virus_coverage_percents = []
let virus_bg_colors = []

let count = 0
Array.prototype.forEach.call(virus_coverage_json, function (virus) {

    if (count < 10)
    {
        virus_coverage_labels.push(virus.virus_name)
        virus_coverage_percents.push(virus.percentage_of_covered)
        virus_bg_colors.push(RandomColor())
    }
    count += 1
})

//virus_coverage_labels.push("Others")
//virus_coverage_percents.push(sum)
//virus_bg_colors.push('rgb'+'(15,10,13)')

const ctx = document.getElementById("virusChart").getContext("2d")

const data = {
  labels: virus_coverage_labels,
  datasets: [{
    label: 'Virus coverage',
    data: virus_coverage_percents,
    backgroundColor: 'rgb(0, 255, 0)',
    hoverOffset: 4
  }]
};

const config = {
  type: 'bar',
  data: data,
  options: {
      plugins: {
          tooltip: {
            callbacks: {
                label: function (context) {
                    console.log(context);
                    let value = context.formattedValue;
                    return value + " %";
                }
            }
          },
      }
  }
};
const virusChart = new Chart(ctx, config);


// --------------------------------------------------------||||||||||||||----------------------------------------------//
// CHarts for all viruses

a = []
Array.prototype.forEach.call(genes_coverage_json, function (virus_with_genes){
    let virus_id = virus_with_genes.virus_id
    let virus_name = virus_with_genes.virus_name
    const ctx = document.getElementById("virus_"+virus_id).getContext("2d")

    let virus_total_length = virus_coverage_json.find( element => element.virus_id === virus_id).from
    //console.log(virus_total_length)
    let gene_coverage_labels = []
    let gene_coverage_percents = []
    let gene_bg_colors = []

    let count = 0
    let sum = 0

    a = virus_with_genes.genes
    a.sort( function (first, second) {
        let first_value = first.from * (first.percentage_of_covered / 100) / virus_total_length * 100
        let second_value = second.from * (second.percentage_of_covered / 100) / virus_total_length * 100
        if (first_value > second_value)
            return -1
        else if (first_value < second_value)
            return 1
        else
            return 0

    });
    Array.prototype.forEach.call(virus_with_genes.genes, function (gene) {
        if (count < 10000 || true)
        {
            let percentage_of_covered = gene.percentage_of_covered
            let size_of_gene = gene.from
            let num_of_covered_of_gene = size_of_gene * (percentage_of_covered / 100)
            let percent_which_covers_gene_in_virus = num_of_covered_of_gene / virus_total_length * 100

            gene_coverage_percents.push(percent_which_covers_gene_in_virus)
            gene_bg_colors.push(RandomColor())
            gene_coverage_labels.push(gene.gene_id+" | "+gene.protein_id)
            sum += percent_which_covers_gene_in_virus
        }
        count += 1
    })

    gene_coverage_labels.push("Not covered")
    gene_coverage_percents.push(100 - sum)
    gene_bg_colors.push('rgb'+'(15,10,13)')

    const data = {
      labels: gene_coverage_labels,
      datasets: [{
        label: 'Virus coverage with genes',
        data: gene_coverage_percents,
        backgroundColor: gene_bg_colors,
        hoverOffset: 4
      }]
    };

    const config = {
      type: 'doughnut',
      data: data,
      options: {
        responsive: true,
        plugins: {
          tooltip: {
            callbacks: {
                label: function (context) {
                    console.log(context);
                    let label = context.label;
                    let value = context.formattedValue;
                    return label + ": " + value + " %";
                }
            }
          },
          legend: {
            position: 'top',
          },
          title: {
            display: true,
            text: virus_name
                }
            }
          },
    };
    console.log("Test");
    const virusWithGenesChart = new Chart(ctx, config);

})





function showMore(id){
    let button = document.getElementById(id)
    let moreRecords = document.getElementsByClassName(id)

    Array.prototype.forEach.call(moreRecords, function (record){
        if (record.style.display === "none")
            record.style.display = "table-row"
        else
            record.style.display = "none"

    })
    if (button.innerHTML === "Show all")
        button.innerHTML = "Show less"
    else
        button.innerHTML = "Show all"
}