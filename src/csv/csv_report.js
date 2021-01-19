$(document).ready(function() {
    $('.table-container').show();
    $('.loading').hide();
    $('#table').bootstrapTable( 'resetView' , {height: window.innerHeight - 200} );
    $('.modal').on('shown.bs.modal', function () {
        window.dispatchEvent(new Event('resize'));
    });
    let to_be_highlighted = window.location.href.toString().split("highlight=").pop();
    let rows = $('.wor');
    rows.each(function() {
        if (this.dataset.index === to_be_highlighted) {
            $(this).children().addClass('active-row');
        }
    });
});
