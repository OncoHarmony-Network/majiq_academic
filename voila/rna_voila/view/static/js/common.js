

$('.close-warnings').click(function(){
    send_ajax('/dismiss-warnings', {})
    .then(() => {
        $('.warnings-container').remove();
    });
});