window.addEventListener('DOMContentLoaded', (event) => {
    main()
});

function main(){
    app = new Vue({
        el: '#app',
        data () {
            return {
                message: "vue loaded.",
                table_info: "table loaded.",
                gff_info: "gff list loaded."
            }
        },
        mounted () {
            axios
                .get('/json')
                .then(response => {this.message = response.data})
            axios
                .get('/table_info')
                .then(response => {this.table_info = response.data})
            axios
                .get('/gff_info')
                .then(response => {this.gff_info = response.data})
        },
        methods : {
            table_relative_path: function(table_path) {
                return table_path.replace(this.table_info.tables_dir_path, '').replace(/^\//, '').replace(/\/$/, '')
            },
            gff_relative_path: function(gff_path) {
                return gff_path.replace(this.gff_info.gff_dir_path, '').replace(/^\//, '').replace(/\/$/, '')
            }
        }
    })
}

function run_pickup(){
    var table_path = document.main_form.table_path.value
    var kmer = document.main_form.kmer.value
    var seqid = document.main_form.seqid.value
    document.pickup_form.table_path.value = table_path
    document.pickup_form.kmer.value = kmer
    document.pickup_form.seqid.value = seqid
    document.pickup_form.submit()
}

function run_lookup(){
    var table_path = document.main_form.table_path.value
    var kmer = document.main_form.kmer.value
    var gff_path = document.main_form.gff_path.value
    var feature_type = document.main_form.feature_type.value
    document.lookup_form.table_path.value = table_path
    document.lookup_form.kmer.value = kmer
    document.lookup_form.gff_path.value = gff_path
    document.lookup_form.feature_type.value = feature_type
    document.lookup_form.submit()
}