$(document).ready(function() {
    $('.table-container').show();
    $('.loading').hide();
    $('#table').bootstrapTable( 'resetView' , {height: window.innerHeight - 200} );
    $('.modal').on('shown.bs.modal', function () {
        window.dispatchEvent(new Event('resize'));
    });
    var i = 0;
    for (const r of data) {
        var table_row = "<tr class=\"wor\" data-idx=\"" + i + "\">";
        i++;
        for (const el of r) {
            var cell = "<td>" + el + "</td>";
            table_row = table_row.concat(cell);
        }
        table_row = table_row.concat("</tr>")
        $(table).find('tbody').append(table_row);
    }
    // Remove "No matching records found" row
    $('table tr.no-records-found').remove();
    $('.fixed-table-border').css("height", "0px");
    let to_be_highlighted = window.location.href.toString().split("highlight=").pop();
    let rows = $('.wor');
    rows.each(function() {
        if (this.dataset.idx === to_be_highlighted) {
            $(this).children().addClass('active-row');
        }
    });
});
