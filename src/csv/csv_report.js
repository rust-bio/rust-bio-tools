$(document).ready(function() {
    $('.table-container').show();
    $('.loading').hide();
    $('#table').bootstrapTable( 'resetView' , {height: window.innerHeight - 200} );
    $('.modal').on('shown.bs.modal', function () {
        window.dispatchEvent(new Event('resize'));
    });
    var decompressed = JSON.parse(LZString.decompressFromUTF16(data));
    let titles = $('#table').bootstrapTable('getVisibleColumns');
    let columns = [];
    for (var x of titles) {
        let title = x.title.split("<a")[0].trim();
        columns.push(title);
    }
    var table_rows = [];
    for (const r of decompressed) {
        var i = 0;
        row = {};
        for (const el of r) {
            row[columns[i]] = el;
            i++;
        }
        table_rows.push(row);
    }
    $('#table').bootstrapTable('append', table_rows)
    let to_be_highlighted = window.location.href.toString().split("highlight=").pop();
    let rows = $("table > tbody > tr");
    rows.each(function() {
        if (this.dataset.index === to_be_highlighted) {
            $(this).children().addClass('active-row');
        }
    });
});
